#include "GeoImage/GeoImage.hpp"
#include "Geographic/GeoHelpers.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "=== Rotated GeoTIFF Test (Central Park) ===\n\n";

    // Central Park center coordinates
    const double centerLat = 40.785;
    const double centerLon = -73.968;
    
    // Image parameters
    const uint32_t imageSize = 512;
    const double imageCenterPixel = imageSize / 2.0;
    
    // Ground coverage: 200x200 meters
    const double groundSizeMeters = 200.0;
    
    // Rotation: 45 degrees
    const double rotationDeg = 45.0;
    const double rotationRad = geohelpers::degreesToRadians(rotationDeg);
    
    // Calculate meters per degree at Central Park latitude
    const double metersPerDegreeLat = geohelpers::EARTH_CIRCUMFERENCE_METERS / 360.0;
    const double metersPerDegreeLon = metersPerDegreeLat * std::cos(geohelpers::degreesToRadians(centerLat));
    
    // Meters per pixel
    const double metersPerPixel = groundSizeMeters / imageSize;
    
    // Degrees per pixel (at image center)
    const double degPerPixelLat = metersPerPixel / metersPerDegreeLat;
    const double degPerPixelLon = metersPerPixel / metersPerDegreeLon;
    
    std::cout << "Image size: " << imageSize << "x" << imageSize << " pixels\n";
    std::cout << "Ground coverage: " << groundSizeMeters << "x" << groundSizeMeters << " meters\n";
    std::cout << "Resolution: " << metersPerPixel << " m/pixel\n";
    std::cout << "Rotation: " << rotationDeg << " degrees\n";
    std::cout << "Center: (" << centerLon << ", " << centerLat << ")\n\n";
    
    // Generate 4 control points at corners of the image
    // We'll compute where each corner pixel maps to in geo coordinates
    // considering the rotation around the center
    
    std::array<double, 4> px, py;  // Pixel coordinates
    std::array<double, 4> gx, gy;  // Geo coordinates (lon, lat)
    
    // Corner pixels (relative to center, then absolute)
    double corners[4][2] = {
        {0, 0},                           // Top-left
        {(double)imageSize, 0},           // Top-right
        {(double)imageSize, (double)imageSize}, // Bottom-right
        {0, (double)imageSize}            // Bottom-left
    };
    
    std::cout << "Control points (pixel -> geo):\n";
    for (int i = 0; i < 4; ++i) {
        px[i] = corners[i][0];
        py[i] = corners[i][1];
        
        // Offset from center pixel
        double dx = px[i] - imageCenterPixel;
        double dy = py[i] - imageCenterPixel;
        
        // Apply rotation
        double rotatedDx = dx * std::cos(rotationRad) - dy * std::sin(rotationRad);
        double rotatedDy = dx * std::sin(rotationRad) + dy * std::cos(rotationRad);
        
        // Convert to geo offset (note: Y is inverted in image coords vs lat)
        double dLon = rotatedDx * degPerPixelLon;
        double dLat = -rotatedDy * degPerPixelLat;  // Negative because pixel Y increases downward
        
        // Final geo coordinates
        gx[i] = centerLon + dLon;
        gy[i] = centerLat + dLat;
        
        std::cout << "  Point " << i << ": (" << px[i] << ", " << py[i] << ") -> ("
                  << gx[i] << ", " << gy[i] << ")\n";
    }
    
    // Create the GeoImage
    geoimage::GeoImage img;
    img.setDimensions(imageSize, imageSize, 3);  // RGB
    
    // Compute and set transformation using the new member function
    if (!img.setTransformationFromPoints(px, py, gx, gy)) {
        std::cout << "\nFAILED to compute affine transformation!\n";
        return 1;
    }
    
    std::cout << "\n4x4 Transformation Matrix:\n";
    const auto& matrix = img.transformation();
    for (int i = 0; i < 4; ++i) {
        std::cout << "  [";
        for (int j = 0; j < 4; ++j) {
            std::cout << std::fixed << std::setprecision(10) << matrix[i * 4 + j];
            if (j < 3) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    
    // Generate a pattern - concentric circles from center
    auto& data = img.data();
    for (uint32_t y = 0; y < imageSize; ++y) {
        for (uint32_t x = 0; x < imageSize; ++x) {
            double dx = x - imageCenterPixel;
            double dy = y - imageCenterPixel;
            double dist = std::sqrt(dx * dx + dy * dy);
            
            // Concentric rings
            float ring = static_cast<float>(std::sin(dist * 0.1) * 0.5 + 0.5);
            
            // Add cross pattern to show rotation
            bool onCross = (std::abs(dx) < 5 || std::abs(dy) < 5);
            
            // Central marker
            bool inCenter = (dist < 10);
            
            size_t idx = (y * imageSize + x) * 3;
            
            if (inCenter) {
                // Red center marker
                data[idx + 0] = 1.0f;
                data[idx + 1] = 0.0f;
                data[idx + 2] = 0.0f;
            } else if (onCross) {
                // Green cross to show orientation
                data[idx + 0] = 0.0f;
                data[idx + 1] = 0.8f;
                data[idx + 2] = 0.0f;
            } else {
                // Blue gradient rings
                data[idx + 0] = ring * 0.3f;
                data[idx + 1] = ring * 0.3f;
                data[idx + 2] = ring * 0.8f;
            }
        }
    }
    
    // Save the GeoTIFF
    std::string filename = "central_park_rotated.tif";
    if (img.save(filename, geoimage::SampleFormat::UInt8)) {
        std::cout << "\nSaved: " << filename << "\n";
    } else {
        std::cout << "\nFAILED to save GeoTIFF!\n";
        return 1;
    }
    
    // Verify by reading back and testing coordinate conversion
    std::cout << "\nVerification:\n";
    geoimage::GeoImage imgRead;
    if (imgRead.open(filename)) {
        // Test center pixel
        double testGeoX, testGeoY;
        imgRead.pixelToGeo(imageCenterPixel, imageCenterPixel, testGeoX, testGeoY);
        std::cout << "  Center pixel (" << imageCenterPixel << ", " << imageCenterPixel 
                  << ") -> (" << testGeoX << ", " << testGeoY << ")\n";
        std::cout << "  Expected: (" << centerLon << ", " << centerLat << ")\n";
        
        double errorLon = std::abs(testGeoX - centerLon);
        double errorLat = std::abs(testGeoY - centerLat);
        std::cout << "  Error: (" << errorLon << ", " << errorLat << ") degrees\n";
        
        // Test resolution
        double mppX, mppY;
        imgRead.getResolutionMetersPerPixel(mppX, mppY);
        std::cout << "  Resolution: " << mppX << " x " << mppY << " m/pixel\n";
        std::cout << "  Expected: ~" << metersPerPixel << " m/pixel\n";
    }
    
    std::cout << "\n=== Test Complete ===\n";
    return 0;
}
