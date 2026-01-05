#include "GeoImage/GeoImage.hpp"
#include <tiffio.h>
#include <algorithm>
#include <cmath>

namespace geoimage {

GeoImage::~GeoImage() {
    close();
}

void GeoImage::setDimensions(uint32_t width, uint32_t height, uint16_t samplesPerPixel) {
    close();
    m_width = width;
    m_height = height;
    m_samplesPerPixel = samplesPerPixel;
    m_bitsPerSample = 32;
    m_data.resize(static_cast<size_t>(width) * height * samplesPerPixel);
    std::fill(m_data.begin(), m_data.end(), 0.0f);
}

GeoImage::GeoImage(GeoImage&& other) noexcept
    : m_tiff(other.m_tiff)
    , m_width(other.m_width)
    , m_height(other.m_height)
    , m_samplesPerPixel(other.m_samplesPerPixel)
    , m_bitsPerSample(other.m_bitsPerSample)
    , m_data(std::move(other.m_data))
    , m_transformation(other.m_transformation)
    , m_hasTransformation(other.m_hasTransformation)
{
    other.m_tiff = nullptr;
    other.m_width = 0;
    other.m_height = 0;
    other.m_hasTransformation = false;
}

GeoImage& GeoImage::operator=(GeoImage&& other) noexcept {
    if (this != &other) {
        close();
        m_tiff = other.m_tiff;
        m_width = other.m_width;
        m_height = other.m_height;
        m_samplesPerPixel = other.m_samplesPerPixel;
        m_bitsPerSample = other.m_bitsPerSample;
        m_data = std::move(other.m_data);
        m_transformation = other.m_transformation;
        m_hasTransformation = other.m_hasTransformation;

        other.m_tiff = nullptr;
        other.m_width = 0;
        other.m_height = 0;
        other.m_hasTransformation = false;
    }
    return *this;
}

bool GeoImage::open(const std::string& filename) {
    close();

    TIFF* tiff = TIFFOpen(filename.c_str(), "r");
    if (!tiff) {
        return false;
    }

    m_tiff = tiff;

    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &m_width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &m_height);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &m_samplesPerPixel);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &m_bitsPerSample);

    // Read GeoTIFF transformation
    m_hasTransformation = false;
    double* transformData = nullptr;
    uint16_t transformCount = 0;
    
    // Try to read ModelTransformationTag (34264)
    if (TIFFGetField(tiff, 34264, &transformCount, &transformData) && transformCount >= 16) {
        for (int i = 0; i < 16; ++i) {
            m_transformation[i] = transformData[i];
        }
        m_hasTransformation = true;
    } else {
        // Try to construct from ModelTiepointTag (33922) and ModelPixelScaleTag (33550)
        double* tiePoints = nullptr;
        double* pixelScale = nullptr;
        uint16_t tieCount = 0, scaleCount = 0;
        
        if (TIFFGetField(tiff, 33922, &tieCount, &tiePoints) && tieCount >= 6 &&
            TIFFGetField(tiff, 33550, &scaleCount, &pixelScale) && scaleCount >= 2) {
            // Construct transformation from tie point and scale
            // TiePoint: [I, J, K, X, Y, Z]
            // PixelScale: [ScaleX, ScaleY, ScaleZ]
            double tpI = tiePoints[0], tpJ = tiePoints[1];
            double tpX = tiePoints[3], tpY = tiePoints[4];
            double scaleX = pixelScale[0];
            double scaleY = pixelScale[1];
            
            setTransformationFromTiePointScale(tpI, tpJ, tpX, tpY, scaleX, scaleY);
        }
    }

    // Read image data as 32-bit float
    m_data.resize(m_width * m_height * m_samplesPerPixel);

    std::vector<uint8_t> rowBuffer(TIFFScanlineSize(tiff));
    for (uint32_t row = 0; row < m_height; ++row) {
        TIFFReadScanline(tiff, rowBuffer.data(), row);
        
        // Convert to float based on source bit depth
        size_t pixelOffset = row * m_width * m_samplesPerPixel;
        for (uint32_t col = 0; col < m_width * m_samplesPerPixel; ++col) {
            if (m_bitsPerSample == 8) {
                m_data[pixelOffset + col] = rowBuffer[col] / 255.0f;
            } else if (m_bitsPerSample == 16) {
                uint16_t val = reinterpret_cast<uint16_t*>(rowBuffer.data())[col];
                m_data[pixelOffset + col] = val / 65535.0f;
            } else if (m_bitsPerSample == 32) {
                m_data[pixelOffset + col] = reinterpret_cast<float*>(rowBuffer.data())[col];
            }
        }
    }

    return true;
}

bool GeoImage::save(const std::string& filename, SampleFormat format) {
    TIFF* tiff = TIFFOpen(filename.c_str(), "w");
    if (!tiff) {
        return false;
    }

    uint16_t bitsPerSample = 8;
    uint16_t sampleFormat = SAMPLEFORMAT_UINT;
    
    switch (format) {
        case SampleFormat::UInt8:
            bitsPerSample = 8;
            sampleFormat = SAMPLEFORMAT_UINT;
            break;
        case SampleFormat::UInt16:
            bitsPerSample = 16;
            sampleFormat = SAMPLEFORMAT_UINT;
            break;
        case SampleFormat::UInt32:
            bitsPerSample = 32;
            sampleFormat = SAMPLEFORMAT_UINT;
            break;
        case SampleFormat::Float32:
            bitsPerSample = 32;
            sampleFormat = SAMPLEFORMAT_IEEEFP;
            break;
    }

    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, m_width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, m_height);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, m_samplesPerPixel);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, bitsPerSample);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, sampleFormat);
    TIFFSetField(tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, m_samplesPerPixel >= 3 ? PHOTOMETRIC_RGB : PHOTOMETRIC_MINISBLACK);

    // Write GeoTIFF ModelTransformationTag (34264) if we have transformation data
    if (m_hasTransformation) {
        TIFFSetField(tiff, 34264, 16, m_transformation.data());
    }

    size_t samplesPerRow = m_width * m_samplesPerPixel;
    
    for (uint32_t row = 0; row < m_height; ++row) {
        const float* rowData = m_data.data() + row * samplesPerRow;
        
        switch (format) {
            case SampleFormat::UInt8: {
                std::vector<uint8_t> buffer(samplesPerRow);
                for (size_t i = 0; i < samplesPerRow; ++i) {
                    float val = std::max(0.0f, std::min(1.0f, rowData[i]));
                    buffer[i] = static_cast<uint8_t>(val * 255.0f);
                }
                TIFFWriteScanline(tiff, buffer.data(), row);
                break;
            }
            case SampleFormat::UInt16: {
                std::vector<uint16_t> buffer(samplesPerRow);
                for (size_t i = 0; i < samplesPerRow; ++i) {
                    float val = std::max(0.0f, std::min(1.0f, rowData[i]));
                    buffer[i] = static_cast<uint16_t>(val * 65535.0f);
                }
                TIFFWriteScanline(tiff, buffer.data(), row);
                break;
            }
            case SampleFormat::UInt32: {
                std::vector<uint32_t> buffer(samplesPerRow);
                for (size_t i = 0; i < samplesPerRow; ++i) {
                    float val = std::max(0.0f, std::min(1.0f, rowData[i]));
                    buffer[i] = static_cast<uint32_t>(val * 4294967295.0);
                }
                TIFFWriteScanline(tiff, buffer.data(), row);
                break;
            }
            case SampleFormat::Float32: {
                TIFFWriteScanline(tiff, const_cast<float*>(rowData), row);
                break;
            }
        }
    }

    TIFFClose(tiff);
    return true;
}

void GeoImage::close() {
    if (m_tiff) {
        TIFFClose(static_cast<TIFF*>(m_tiff));
        m_tiff = nullptr;
    }
}

void GeoImage::clear() {
    close();
    m_width = 0;
    m_height = 0;
    m_samplesPerPixel = 1;
    m_bitsPerSample = 32;
    m_data.clear();
    m_transformation = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    m_hasTransformation = false;
}

void GeoImage::setTransformation(const ModelTransformation& transform) {
    m_transformation = transform;
    m_hasTransformation = true;
}

void GeoImage::setTransformationFromTiePointScale(double tiePointX, double tiePointY,
                                                   double geoX, double geoY,
                                                   double scaleX, double scaleY) {
    // Build 4x4 transformation matrix from tie point and pixel scale
    // X = geoX + (pixelX - tiePointX) * scaleX
    // Y = geoY - (pixelY - tiePointY) * scaleY  (Y typically inverted)
    //
    // Matrix form (row-major):
    // | scaleX    0      0    geoX - tiePointX * scaleX  |
    // |   0    -scaleY   0    geoY + tiePointY * scaleY  |
    // |   0       0      0              0                |
    // |   0       0      0              1                |
    
    m_transformation = {
        scaleX,   0.0,      0.0,  geoX - tiePointX * scaleX,
        0.0,     -scaleY,   0.0,  geoY + tiePointY * scaleY,
        0.0,      0.0,      0.0,  0.0,
        0.0,      0.0,      0.0,  1.0
    };
    m_hasTransformation = true;
}

void GeoImage::pixelToGeo(double pixelX, double pixelY, double& geoX, double& geoY) const {
    // Apply transformation matrix: [x, y, z, 1] = M * [i, j, k, 1]
    // For 2D: x = m[0]*i + m[1]*j + m[3]
    //         y = m[4]*i + m[5]*j + m[7]
    geoX = m_transformation[0] * pixelX + m_transformation[1] * pixelY + m_transformation[3];
    geoY = m_transformation[4] * pixelX + m_transformation[5] * pixelY + m_transformation[7];
}

void GeoImage::geoToPixel(double geoX, double geoY, double& pixelX, double& pixelY) const {
    // Inverse of 2D affine transformation
    // For simple scale+translate: pixelX = (geoX - tx) / scaleX
    //                             pixelY = (geoY - ty) / scaleY
    double a = m_transformation[0];  // scaleX
    double b = m_transformation[1];
    double e = m_transformation[4];
    double f = m_transformation[5];  // -scaleY
    double tx = m_transformation[3];
    double ty = m_transformation[7];
    
    double det = a * f - b * e;
    if (std::abs(det) > 1e-10) {
        pixelX = (f * (geoX - tx) - b * (geoY - ty)) / det;
        pixelY = (-e * (geoX - tx) + a * (geoY - ty)) / det;
    } else {
        pixelX = 0.0;
        pixelY = 0.0;
    }
}

} // namespace geoimage
