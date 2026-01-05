#include "GeoImage/GeoImage.hpp"
#include "Geographic/GeoHelpers.hpp"
#include "BLAS/blas.hpp"
#include <tiffio.h>
#include <algorithm>
#include <cmath>

namespace geoimage {

// GeoTIFF tag definitions for libtiff
static const TIFFFieldInfo geotiffFieldInfo[] = {
    { 33550, -1, -1, TIFF_DOUBLE, FIELD_CUSTOM, true, true, const_cast<char*>("ModelPixelScaleTag") },
    { 33922, -1, -1, TIFF_DOUBLE, FIELD_CUSTOM, true, true, const_cast<char*>("ModelTiepointTag") },
    { 34264, -1, -1, TIFF_DOUBLE, FIELD_CUSTOM, true, true, const_cast<char*>("ModelTransformationTag") },
    { 34735, -1, -1, TIFF_SHORT, FIELD_CUSTOM, true, true, const_cast<char*>("GeoKeyDirectoryTag") },
    { 34736, -1, -1, TIFF_DOUBLE, FIELD_CUSTOM, true, true, const_cast<char*>("GeoDoubleParamsTag") },
    { 34737, -1, -1, TIFF_ASCII, FIELD_CUSTOM, true, false, const_cast<char*>("GeoAsciiParamsTag") },
};

static void registerGeoTIFFTags(TIFF* tiff) {
    TIFFMergeFieldInfo(tiff, geotiffFieldInfo, sizeof(geotiffFieldInfo) / sizeof(geotiffFieldInfo[0]));
}

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

    // Register GeoTIFF tags
    registerGeoTIFFTags(tiff);

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
    if (TIFFGetField(tiff, geohelpers::TIFFTAG_MODELTRANSFORMATION, &transformCount, &transformData) && transformCount >= 16) {
        for (int i = 0; i < 16; ++i) {
            m_transformation[i] = transformData[i];
        }
        m_hasTransformation = true;
    } else {
        // Try to construct from ModelTiepointTag (33922) and ModelPixelScaleTag (33550)
        double* tiePoints = nullptr;
        double* pixelScale = nullptr;
        uint16_t tieCount = 0, scaleCount = 0;
        
        if (TIFFGetField(tiff, geohelpers::TIFFTAG_MODELTIEPOINT, &tieCount, &tiePoints) && tieCount >= 6 &&
            TIFFGetField(tiff, geohelpers::TIFFTAG_MODELPIXELSCALE, &scaleCount, &pixelScale) && scaleCount >= 2) {
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

    // Register GeoTIFF tags
    registerGeoTIFFTags(tiff);

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

    // Write GeoTIFF tags if we have transformation data
    if (m_hasTransformation) {
        // ModelTransformationTag (34264)
        TIFFSetField(tiff, geohelpers::TIFFTAG_MODELTRANSFORMATION, 16, m_transformation.data());
        
        // GeoKeyDirectoryTag (34735) for EPSG:4326
        TIFFSetField(tiff, geohelpers::TIFFTAG_GEOKEYDIRECTORY, 
                     geohelpers::EPSG_4326_GeoKeyDirectory.size(), 
                     geohelpers::EPSG_4326_GeoKeyDirectory.data());
        
        // GeoDoubleParamsTag (34736) - WGS84 ellipsoid parameters
        TIFFSetField(tiff, geohelpers::TIFFTAG_GEODOUBLEPARAMS, 
                     geohelpers::EPSG_4326_GeoDoubleParams.size(), 
                     geohelpers::EPSG_4326_GeoDoubleParams.data());
        
        // Calculate and set XResolution/YResolution in pixels per centimeter
        double metersPerPixelX, metersPerPixelY;
        getResolutionMetersPerPixel(metersPerPixelX, metersPerPixelY);
        if (metersPerPixelX > 0 && metersPerPixelY > 0) {
            // Convert meters per pixel to pixels per centimeter
            // pixels/cm = 0.01 / meters_per_pixel
            float xRes = static_cast<float>(0.01 / metersPerPixelX);
            float yRes = static_cast<float>(0.01 / metersPerPixelY);
            TIFFSetField(tiff, TIFFTAG_XRESOLUTION, xRes);
            TIFFSetField(tiff, TIFFTAG_YRESOLUTION, yRes);
            TIFFSetField(tiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_CENTIMETER);
        }
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

void GeoImage::getResolutionMetersPerPixel(double& metersPerPixelX, double& metersPerPixelY) const {
    if (!m_hasTransformation) {
        metersPerPixelX = 0.0;
        metersPerPixelY = 0.0;
        return;
    }
    
    // Get the center of the image in geo coordinates
    double centerPixelX = m_width / 2.0;
    double centerPixelY = m_height / 2.0;
    double centerGeoX, centerGeoY;
    pixelToGeo(centerPixelX, centerPixelY, centerGeoX, centerGeoY);
    
    // Get pixel scale from transformation matrix (degrees per pixel)
    double degPerPixelX = std::abs(m_transformation[0]);  // scaleX
    double degPerPixelY = std::abs(m_transformation[5]);  // scaleY (may be negative)
    
    // Convert degrees to meters at the center latitude
    // 1 degree of latitude ≈ 111,320 meters (constant)
    // 1 degree of longitude ≈ 111,320 * cos(latitude) meters
    double latRad = geohelpers::degreesToRadians(centerGeoY);
    double metersPerDegreeLat = geohelpers::EARTH_CIRCUMFERENCE_METERS / 360.0;
    double metersPerDegreeLon = metersPerDegreeLat * std::cos(latRad);
    
    metersPerPixelX = degPerPixelX * metersPerDegreeLon;
    metersPerPixelY = degPerPixelY * metersPerDegreeLat;
}

bool GeoImage::setTransformationFromPoints(const std::array<double, 3>& pixelX, const std::array<double, 3>& pixelY,
                                           const std::array<double, 3>& geoX, const std::array<double, 3>& geoY) {
    blas::Vector<6> coeffs;
    if (!blas::affineFrom3Points(pixelX, pixelY, geoX, geoY, coeffs)) {
        return false;
    }
    auto matrix = blas::affineCoeffsToMatrix4x4(coeffs);
    for (int i = 0; i < 16; ++i) {
        m_transformation[i] = matrix[i];
    }
    m_hasTransformation = true;
    return true;
}

bool GeoImage::setTransformationFromPoints(const std::array<double, 4>& pixelX, const std::array<double, 4>& pixelY,
                                           const std::array<double, 4>& geoX, const std::array<double, 4>& geoY) {
    blas::Vector<6> coeffs;
    if (!blas::affineFrom4Points(pixelX, pixelY, geoX, geoY, coeffs)) {
        return false;
    }
    auto matrix = blas::affineCoeffsToMatrix4x4(coeffs);
    for (int i = 0; i < 16; ++i) {
        m_transformation[i] = matrix[i];
    }
    m_hasTransformation = true;
    return true;
}

bool GeoImage::setTransformationFromPoints(const std::array<double, 5>& pixelX, const std::array<double, 5>& pixelY,
                                           const std::array<double, 5>& geoX, const std::array<double, 5>& geoY) {
    blas::Vector<6> coeffs;
    if (!blas::affineFrom5Points(pixelX, pixelY, geoX, geoY, coeffs)) {
        return false;
    }
    auto matrix = blas::affineCoeffsToMatrix4x4(coeffs);
    for (int i = 0; i < 16; ++i) {
        m_transformation[i] = matrix[i];
    }
    m_hasTransformation = true;
    return true;
}

} // namespace geoimage
