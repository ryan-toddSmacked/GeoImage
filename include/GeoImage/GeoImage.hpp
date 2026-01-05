#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace geoimage {

class GeoImage {
public:
    GeoImage() = default;
    ~GeoImage();

    // Disable copy, allow move
    GeoImage(const GeoImage&) = delete;
    GeoImage& operator=(const GeoImage&) = delete;
    GeoImage(GeoImage&&) noexcept;
    GeoImage& operator=(GeoImage&&) noexcept;

    // File operations
    bool open(const std::string& filename);
    bool save(const std::string& filename);
    void close();

    // Image properties
    uint32_t width() const { return m_width; }
    uint32_t height() const { return m_height; }
    uint16_t samplesPerPixel() const { return m_samplesPerPixel; }
    uint16_t bitsPerSample() const { return m_bitsPerSample; }

    // Data access
    const std::vector<uint8_t>& data() const { return m_data; }
    std::vector<uint8_t>& data() { return m_data; }

    bool isOpen() const { return m_tiff != nullptr; }

private:
    void* m_tiff = nullptr;  // TIFF* handle
    uint32_t m_width = 0;
    uint32_t m_height = 0;
    uint16_t m_samplesPerPixel = 1;
    uint16_t m_bitsPerSample = 8;
    std::vector<uint8_t> m_data;
};

} // namespace geoimage
