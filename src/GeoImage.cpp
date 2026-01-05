#include "GeoImage/GeoImage.hpp"
#include <tiffio.h>

namespace geoimage {

GeoImage::~GeoImage() {
    close();
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

    // Read image data
    size_t bytesPerPixel = (m_bitsPerSample / 8) * m_samplesPerPixel;
    m_data.resize(m_width * m_height * bytesPerPixel);

    for (uint32_t row = 0; row < m_height; ++row) {
        TIFFReadScanline(tiff, m_data.data() + row * m_width * bytesPerPixel, row);
    }

    return true;
}

bool GeoImage::save(const std::string& filename) {
    TIFF* tiff = TIFFOpen(filename.c_str(), "w");
    if (!tiff) {
        return false;
    }

    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, m_width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, m_height);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, m_samplesPerPixel);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, m_bitsPerSample);
    TIFFSetField(tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

    size_t bytesPerPixel = (m_bitsPerSample / 8) * m_samplesPerPixel;
    for (uint32_t row = 0; row < m_height; ++row) {
        TIFFWriteScanline(tiff, m_data.data() + row * m_width * bytesPerPixel, row);
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
