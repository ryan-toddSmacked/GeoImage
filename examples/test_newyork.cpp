#include "GeoImage/GeoImage.hpp"
#include "Geographic/GeoHelpers.hpp"
#include <iostream>
#include <cmath>

int main() {
    std::cout << "=== New York GeoTIFF Test ===\n\n";

    // New York City area coordinates
    // Approximate bounding box around Manhattan
    const double minLon = -74.05;  // West
    const double maxLon = -73.90;  // East
    const double minLat = 40.70;   // South
    const double maxLat = 40.85;   // North

    // Image dimensions
    const uint32_t width = 512;
    const uint32_t height = 512;
    const uint16_t samples = 3;  // RGB

    // Calculate pixel scale (degrees per pixel)
    double scaleX = (maxLon - minLon) / width;
    double scaleY = (maxLat - minLat) / height;

    std::cout << "Image dimensions: " << width << "x" << height << " pixels\n";
    std::cout << "Geographic extent:\n";
    std::cout << "  Longitude: " << minLon << " to " << maxLon << " (" << (maxLon - minLon) << " degrees)\n";
    std::cout << "  Latitude:  " << minLat << " to " << maxLat << " (" << (maxLat - minLat) << " degrees)\n";
    std::cout << "  Scale X: " << scaleX << " deg/pixel\n";
    std::cout << "  Scale Y: " << scaleY << " deg/pixel\n\n";

    // Create the image
    geoimage::GeoImage img;
    img.setDimensions(width, height, samples);

    // Set the geotransformation
    // Top-left corner is (minLon, maxLat) - GeoTIFF convention
    img.setTransformationFromTiePointScale(0, 0, minLon, maxLat, scaleX, scaleY);

    auto& data = img.data();

    // Generate a pattern that represents terrain-like data
    // Using a combination of gradients and patterns
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {
            // Convert pixel to geo coordinates
            double geoX, geoY;
            img.pixelToGeo(static_cast<double>(x), static_cast<double>(y), geoX, geoY);

            // Create a pattern based on distance from Central Park (approx 40.785, -73.968)
            double centralParkLat = 40.785;
            double centralParkLon = -73.968;
            
            double distLat = geoY - centralParkLat;
            double distLon = geoX - centralParkLon;
            double dist = std::sqrt(distLat * distLat + distLon * distLon);

            // Create concentric rings emanating from Central Park
            float ring = std::sin(dist * 200.0) * 0.5f + 0.5f;
            
            // Add some "water" blue for the Hudson River (west side)
            bool isHudson = (geoX < -74.01);
            
            // Add some "water" blue for the East River
            bool isEastRiver = (geoX > -73.93 && geoY < 40.78);

            size_t idx = (y * width + x) * samples;
            
            if (isHudson || isEastRiver) {
                // Water - blue
                data[idx + 0] = 0.1f;   // R
                data[idx + 1] = 0.3f;   // G
                data[idx + 2] = 0.7f;   // B
            } else if (dist < 0.02) {
                // Central Park area - green
                float green = 0.4f + ring * 0.3f;
                data[idx + 0] = 0.2f;
                data[idx + 1] = green;
                data[idx + 2] = 0.2f;
            } else {
                // Urban area - gray with pattern
                float urban = 0.3f + ring * 0.4f;
                data[idx + 0] = urban;
                data[idx + 1] = urban * 0.9f;
                data[idx + 2] = urban * 0.85f;
            }
        }
    }

    // Save as different formats
    std::string basePath = "newyork_geotiff";
    
    std::cout << "Saving GeoTIFF files...\n";
    
    if (img.save(basePath + "_u8.tif", geoimage::SampleFormat::UInt8)) {
        std::cout << "  Saved: " << basePath << "_u8.tif (8-bit)\n";
    } else {
        std::cout << "  FAILED to save 8-bit version\n";
    }
    
    if (img.save(basePath + "_u16.tif", geoimage::SampleFormat::UInt16)) {
        std::cout << "  Saved: " << basePath << "_u16.tif (16-bit)\n";
    } else {
        std::cout << "  FAILED to save 16-bit version\n";
    }
    
    if (img.save(basePath + "_f32.tif", geoimage::SampleFormat::Float32)) {
        std::cout << "  Saved: " << basePath << "_f32.tif (32-bit float)\n";
    } else {
        std::cout << "  FAILED to save float version\n";
    }

    // Test reading it back
    std::cout << "\nReading back the GeoTIFF...\n";
    geoimage::GeoImage imgRead;
    if (imgRead.open(basePath + "_u8.tif")) {
        std::cout << "  Dimensions: " << imgRead.width() << "x" << imgRead.height() << "\n";
        std::cout << "  Samples per pixel: " << imgRead.samplesPerPixel() << "\n";
        std::cout << "  Has transformation: " << (imgRead.hasTransformation() ? "yes" : "no") << "\n";
        
        if (imgRead.hasTransformation()) {
            // Test coordinate conversion
            double testGeoX, testGeoY;
            imgRead.pixelToGeo(0, 0, testGeoX, testGeoY);
            std::cout << "  Top-left corner: (" << testGeoX << ", " << testGeoY << ")\n";
            
            imgRead.pixelToGeo(width - 1, height - 1, testGeoX, testGeoY);
            std::cout << "  Bottom-right corner: (" << testGeoX << ", " << testGeoY << ")\n";
            
            // Output resolution in meters per pixel
            double metersPerPixelX, metersPerPixelY;
            imgRead.getResolutionMetersPerPixel(metersPerPixelX, metersPerPixelY);
            std::cout << "  Resolution X: " << metersPerPixelX << " m/pixel\n";
            std::cout << "  Resolution Y: " << metersPerPixelY << " m/pixel\n";
            
            // Test reverse conversion
            double testPixelX, testPixelY;
            imgRead.geoToPixel(-73.968, 40.785, testPixelX, testPixelY);
            std::cout << "  Central Park (-73.968, 40.785) -> pixel (" 
                      << testPixelX << ", " << testPixelY << ")\n";
        }
    } else {
        std::cout << "  FAILED to read back GeoTIFF\n";
    }

    std::cout << "\n=== Test Complete ===\n";
    return 0;
}
