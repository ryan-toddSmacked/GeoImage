#include "GeoImage/GeoImage.hpp"
#include <iostream>
#include <cmath>
#include <cassert>
#include <filesystem>

namespace fs = std::filesystem;

// Test output directory
const std::string TEST_OUTPUT_DIR = "test_output";

// Helper to ensure test output directory exists
void ensureOutputDir() {
    fs::create_directories(TEST_OUTPUT_DIR);
}

// Generate a horizontal gradient histogram image
bool testHorizontalGradient() {
    std::cout << "Test: Horizontal Gradient... ";
    
    geoimage::GeoImage img;
    const uint32_t width = 256;
    const uint32_t height = 64;
    
    // Create grayscale gradient
    img.setDimensions(width, height, 1);
    auto& data = img.data();
    
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {
            data[y * width + x] = static_cast<float>(x) / (width - 1);
        }
    }
    
    // Save in different formats
    std::string basePath = TEST_OUTPUT_DIR + "/gradient_h";
    if (!img.save(basePath + "_u8.tif", geoimage::SampleFormat::UInt8)) return false;
    if (!img.save(basePath + "_u16.tif", geoimage::SampleFormat::UInt16)) return false;
    if (!img.save(basePath + "_f32.tif", geoimage::SampleFormat::Float32)) return false;
    
    std::cout << "PASSED\n";
    return true;
}

// Generate a vertical gradient histogram image
bool testVerticalGradient() {
    std::cout << "Test: Vertical Gradient... ";
    
    geoimage::GeoImage img;
    const uint32_t width = 64;
    const uint32_t height = 256;
    
    img.setDimensions(width, height, 1);
    auto& data = img.data();
    
    for (uint32_t y = 0; y < height; ++y) {
        float value = static_cast<float>(y) / (height - 1);
        for (uint32_t x = 0; x < width; ++x) {
            data[y * width + x] = value;
        }
    }
    
    if (!img.save(TEST_OUTPUT_DIR + "/gradient_v_u8.tif", geoimage::SampleFormat::UInt8)) return false;
    
    std::cout << "PASSED\n";
    return true;
}

// Generate RGB color bars
bool testColorBars() {
    std::cout << "Test: RGB Color Bars... ";
    
    geoimage::GeoImage img;
    const uint32_t width = 256;
    const uint32_t height = 64;
    const uint16_t samples = 3; // RGB
    
    img.setDimensions(width, height, samples);
    auto& data = img.data();
    
    // Create 8 color bars: Black, Red, Green, Yellow, Blue, Magenta, Cyan, White
    const float colors[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}
    };
    
    uint32_t barWidth = width / 8;
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {
            int colorIdx = std::min(static_cast<int>(x / barWidth), 7);
            size_t idx = (y * width + x) * samples;
            data[idx + 0] = colors[colorIdx][0];
            data[idx + 1] = colors[colorIdx][1];
            data[idx + 2] = colors[colorIdx][2];
        }
    }
    
    if (!img.save(TEST_OUTPUT_DIR + "/color_bars_u8.tif", geoimage::SampleFormat::UInt8)) return false;
    if (!img.save(TEST_OUTPUT_DIR + "/color_bars_u16.tif", geoimage::SampleFormat::UInt16)) return false;
    
    std::cout << "PASSED\n";
    return true;
}

// Generate a histogram distribution image
bool testHistogramImage() {
    std::cout << "Test: Histogram Distribution... ";
    
    geoimage::GeoImage img;
    const uint32_t width = 256;
    const uint32_t height = 128;
    
    img.setDimensions(width, height, 1);
    auto& data = img.data();
    
    // Create a Gaussian-like histogram visualization
    for (uint32_t x = 0; x < width; ++x) {
        float normalized = static_cast<float>(x) / width;
        float center = 0.5f;
        float sigma = 0.15f;
        float gaussVal = std::exp(-0.5f * std::pow((normalized - center) / sigma, 2.0f));
        uint32_t barHeight = static_cast<uint32_t>(gaussVal * (height - 1));
        
        for (uint32_t y = 0; y < height; ++y) {
            uint32_t flippedY = height - 1 - y;
            data[y * width + x] = (flippedY < barHeight) ? 1.0f : 0.1f;
        }
    }
    
    if (!img.save(TEST_OUTPUT_DIR + "/histogram_gauss_u8.tif", geoimage::SampleFormat::UInt8)) return false;
    
    std::cout << "PASSED\n";
    return true;
}

// Test read/write cycle - save and reopen
bool testReadWriteCycle() {
    std::cout << "Test: Read/Write Cycle... ";
    
    geoimage::GeoImage imgWrite;
    const uint32_t width = 100;
    const uint32_t height = 100;
    
    imgWrite.setDimensions(width, height, 1);
    auto& writeData = imgWrite.data();
    
    // Fill with checkerboard pattern
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {
            bool isWhite = ((x / 10) + (y / 10)) % 2 == 0;
            writeData[y * width + x] = isWhite ? 1.0f : 0.0f;
        }
    }
    
    std::string filepath = TEST_OUTPUT_DIR + "/checkerboard_u8.tif";
    if (!imgWrite.save(filepath, geoimage::SampleFormat::UInt8)) {
        std::cout << "FAILED (save)\n";
        return false;
    }
    
    // Read it back
    geoimage::GeoImage imgRead;
    if (!imgRead.open(filepath)) {
        std::cout << "FAILED (open)\n";
        return false;
    }
    
    // Verify dimensions
    if (imgRead.width() != width || imgRead.height() != height) {
        std::cout << "FAILED (dimensions mismatch)\n";
        return false;
    }
    
    // Verify data (with tolerance for float conversion)
    const auto& readData = imgRead.data();
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {
            float expected = writeData[y * width + x];
            float actual = readData[y * width + x];
            if (std::abs(expected - actual) > 0.01f) {
                std::cout << "FAILED (data mismatch at " << x << "," << y << ")\n";
                return false;
            }
        }
    }
    
    std::cout << "PASSED\n";
    return true;
}

// Test image operations - invert
bool testInvertOperation() {
    std::cout << "Test: Invert Operation... ";
    
    geoimage::GeoImage img;
    const uint32_t width = 128;
    const uint32_t height = 128;
    
    img.setDimensions(width, height, 1);
    auto& data = img.data();
    
    // Create gradient
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {
            data[y * width + x] = static_cast<float>(x) / (width - 1);
        }
    }
    
    // Save original
    img.save(TEST_OUTPUT_DIR + "/invert_before.tif", geoimage::SampleFormat::UInt8);
    
    // Invert the image
    for (auto& val : data) {
        val = 1.0f - val;
    }
    
    // Save inverted
    img.save(TEST_OUTPUT_DIR + "/invert_after.tif", geoimage::SampleFormat::UInt8);
    
    // Verify inversion
    if (std::abs(data[0] - 1.0f) > 0.01f || std::abs(data[width - 1] - 0.0f) > 0.01f) {
        std::cout << "FAILED (inversion incorrect)\n";
        return false;
    }
    
    std::cout << "PASSED\n";
    return true;
}

// Test brightness/contrast adjustment
bool testBrightnessContrast() {
    std::cout << "Test: Brightness/Contrast Adjustment... ";
    
    geoimage::GeoImage img;
    const uint32_t width = 256;
    const uint32_t height = 64;
    
    img.setDimensions(width, height, 1);
    auto& data = img.data();
    
    // Create gradient
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {
            data[y * width + x] = static_cast<float>(x) / (width - 1);
        }
    }
    
    img.save(TEST_OUTPUT_DIR + "/bc_original.tif", geoimage::SampleFormat::UInt8);
    
    // Increase contrast (multiply by 1.5, centered at 0.5)
    for (auto& val : data) {
        val = (val - 0.5f) * 1.5f + 0.5f;
        val = std::max(0.0f, std::min(1.0f, val)); // Clamp
    }
    
    img.save(TEST_OUTPUT_DIR + "/bc_high_contrast.tif", geoimage::SampleFormat::UInt8);
    
    std::cout << "PASSED\n";
    return true;
}

// Test float precision preservation
bool testFloatPrecision() {
    std::cout << "Test: Float32 Precision... ";
    
    geoimage::GeoImage imgWrite;
    const uint32_t width = 10;
    const uint32_t height = 10;
    
    imgWrite.setDimensions(width, height, 1);
    auto& writeData = imgWrite.data();
    
    // Use specific float values that need precision
    float testValues[] = {0.0f, 0.123456f, 0.5f, 0.789012f, 1.0f, 1.5f, -0.5f, 100.0f};
    for (uint32_t i = 0; i < width * height && i < 8; ++i) {
        writeData[i] = testValues[i % 8];
    }
    
    std::string filepath = TEST_OUTPUT_DIR + "/precision_f32.tif";
    if (!imgWrite.save(filepath, geoimage::SampleFormat::Float32)) {
        std::cout << "FAILED (save)\n";
        return false;
    }
    
    geoimage::GeoImage imgRead;
    if (!imgRead.open(filepath)) {
        std::cout << "FAILED (open)\n";
        return false;
    }
    
    const auto& readData = imgRead.data();
    for (uint32_t i = 0; i < 8; ++i) {
        if (std::abs(writeData[i] - readData[i]) > 1e-6f) {
            std::cout << "FAILED (precision loss at index " << i << ")\n";
            return false;
        }
    }
    
    std::cout << "PASSED\n";
    return true;
}

int main() {
    std::cout << "=== GeoImage Test Suite ===\n\n";
    
    ensureOutputDir();
    
    int passed = 0;
    int failed = 0;
    
    auto runTest = [&](bool (*testFunc)()) {
        if (testFunc()) {
            ++passed;
        } else {
            ++failed;
        }
    };
    
    runTest(testHorizontalGradient);
    runTest(testVerticalGradient);
    runTest(testColorBars);
    runTest(testHistogramImage);
    runTest(testReadWriteCycle);
    runTest(testInvertOperation);
    runTest(testBrightnessContrast);
    runTest(testFloatPrecision);
    
    std::cout << "\n=== Results ===\n";
    std::cout << "Passed: " << passed << "\n";
    std::cout << "Failed: " << failed << "\n";
    
    std::cout << "\nTest images saved to: " << fs::absolute(TEST_OUTPUT_DIR) << "\n";
    
    return failed > 0 ? 1 : 0;
}
