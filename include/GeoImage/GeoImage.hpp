#pragma once

#include <string>
#include <vector>
#include <array>
#include <cstdint>

namespace geoimage {

enum class SampleFormat {
    UInt8,
    UInt16,
    UInt32,
    Float32
};

// 4x4 affine transformation matrix for GeoTIFF (stored row-major)
// Transforms raster (i, j, k) to model coordinates (x, y, z)
// | x |   | a  b  c  d |   | i |
// | y | = | e  f  g  h | * | j |
// | z |   | i  j  k  l |   | k |
// | 1 |   | m  n  o  p |   | 1 |
using ModelTransformation = std::array<double, 16>;

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
    void clear();

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

    // GeoTIFF transformation matrix
    const ModelTransformation& transformation() const { return m_transformation; }
    void setTransformation(const ModelTransformation& transform);
    bool hasTransformation() const { return m_hasTransformation; }
    
    // Create transformation from tie point and pixel scale
    void setTransformationFromTiePointScale(double tiePointX, double tiePointY,
                                            double geoX, double geoY,
                                            double scaleX, double scaleY);
    
    // Convert pixel coordinates to geo coordinates using transformation
    void pixelToGeo(double pixelX, double pixelY, double& geoX, double& geoY) const;
    void geoToPixel(double geoX, double geoY, double& pixelX, double& pixelY) const;

    bool isOpen() const { return m_tiff != nullptr; }

private:
    void* m_tiff = nullptr;  // TIFF* handle
    uint32_t m_width = 0;
    uint32_t m_height = 0;
    uint16_t m_samplesPerPixel = 1;
    uint16_t m_bitsPerSample = 32;
    std::vector<float> m_data;
    
    // GeoTIFF data
    ModelTransformation m_transformation = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};  // Identity
    bool m_hasTransformation = false;
};

} // namespace geoimage
