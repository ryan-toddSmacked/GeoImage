#include "GeoImage/GeoImage.hpp"
#include <tiffio.h>
#include <algorithm>

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
{
    other.m_tiff = nullptr;
    other.m_width = 0;
    other.m_height = 0;
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

        other.m_tiff = nullptr;
        other.m_width = 0;
        other.m_height = 0;
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
    m_width = 0;
    m_height = 0;
    m_data.clear();
}

} // namespace geoimage
