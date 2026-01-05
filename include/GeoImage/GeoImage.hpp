#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace geoimage {

enum class SampleFormat {
    UInt8,
    UInt16,
    UInt32,
    Float32
};

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
    bool save(const std::string& filename, SampleFormat format = SampleFormat::UInt8);
    void close();

    // Image properties
    uint32_t width() const { return m_width; }
    uint32_t height() const { return m_height; }
    uint16_t samplesPerPixel() const { return m_samplesPerPixel; }
    uint16_t bitsPerSample() const { return m_bitsPerSample; }

    // Create/resize image
    void setDimensions(uint32_t width, uint32_t height, uint16_t samplesPerPixel = 1);

    // Data access
    const std::vector<float>& data() const { return m_data; }
    std::vector<float>& data() { return m_data; }

    bool isOpen() const { return m_tiff != nullptr; }

private:
    void* m_tiff = nullptr;  // TIFF* handle
    uint32_t m_width = 0;
    uint32_t m_height = 0;
    uint16_t m_samplesPerPixel = 1;
    uint16_t m_bitsPerSample = 32;
    std::vector<float> m_data;
};

} // namespace geoimage
