#pragma once

#include <cmath>
#include <vector>
#include <array>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace geohelpers {

    // WGS84 ellipsoid constants
    constexpr double EARTH_RADIUS_METERS = 6378137.0;              // Semi-major axis (equatorial radius)
    constexpr double EARTH_RADIUS_POLAR_METERS = 6356752.314245;   // Semi-minor axis (polar radius)
    constexpr double EARTH_CIRCUMFERENCE_METERS = 2.0 * M_PI * EARTH_RADIUS_METERS;
    constexpr double EARTH_FLATTENING = 1.0 / 298.257223563;       // WGS84 flattening

    // Unit conversion constants
    constexpr double METERS_PER_FOOT = 0.3048;
    constexpr double FEET_PER_METER = 1.0 / METERS_PER_FOOT;

    //--------------------------------------------------------------------------
    // Angle conversions
    //--------------------------------------------------------------------------

    constexpr double degreesToRadians(double degrees) {
        return degrees * (M_PI / 180.0);
    }

    constexpr double radiansToDegrees(double radians) {
        return radians * (180.0 / M_PI);
    }

    //--------------------------------------------------------------------------
    // Coordinate system conversions (Web Mercator / EPSG:3857)
    //--------------------------------------------------------------------------

    // Convert lat/lon (EPSG:4326) to Web Mercator meters (EPSG:3857)
    inline void latLonToMercator(double lat, double lon, double& outX, double& outY) {
        outX = EARTH_RADIUS_METERS * degreesToRadians(lon);
        outY = EARTH_RADIUS_METERS * std::log(std::tan(M_PI / 4.0 + degreesToRadians(lat) / 2.0));
    }

    // Convert Web Mercator meters (EPSG:3857) to lat/lon (EPSG:4326)
    inline void mercatorToLatLon(double x, double y, double& outLat, double& outLon) {
        outLon = radiansToDegrees(x / EARTH_RADIUS_METERS);
        outLat = radiansToDegrees(2.0 * std::atan(std::exp(y / EARTH_RADIUS_METERS)) - M_PI / 2.0);
    }

    //--------------------------------------------------------------------------
    // Distance calculations
    //--------------------------------------------------------------------------

    // Haversine formula - great circle distance between two lat/lon points
    inline double haversineDistance(double lat1, double lon1, double lat2, double lon2) {
        double dLat = degreesToRadians(lat2 - lat1);
        double dLon = degreesToRadians(lon2 - lon1);
        
        double a = std::sin(dLat / 2.0) * std::sin(dLat / 2.0) +
                   std::cos(degreesToRadians(lat1)) * std::cos(degreesToRadians(lat2)) *
                   std::sin(dLon / 2.0) * std::sin(dLon / 2.0);
        
        double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
        
        return EARTH_RADIUS_METERS * c;
    }

    // Calculate bearing from point 1 to point 2 (in degrees, 0 = North, 90 = East)
    inline double bearing(double lat1, double lon1, double lat2, double lon2) {
        double dLon = degreesToRadians(lon2 - lon1);
        double lat1Rad = degreesToRadians(lat1);
        double lat2Rad = degreesToRadians(lat2);
        
        double y = std::sin(dLon) * std::cos(lat2Rad);
        double x = std::cos(lat1Rad) * std::sin(lat2Rad) -
                   std::sin(lat1Rad) * std::cos(lat2Rad) * std::cos(dLon);
        
        double bearingRad = std::atan2(y, x);
        double bearingDeg = radiansToDegrees(bearingRad);
        
        // Normalize to 0-360
        return std::fmod(bearingDeg + 360.0, 360.0);
    }

    // Calculate destination point given start point, bearing, and distance
    inline void destinationPoint(double lat, double lon, double bearingDeg, double distanceMeters,
                                 double& outLat, double& outLon) {
        double bearingRad = degreesToRadians(bearingDeg);
        double latRad = degreesToRadians(lat);
        double lonRad = degreesToRadians(lon);
        double angularDist = distanceMeters / EARTH_RADIUS_METERS;
        
        double outLatRad = std::asin(std::sin(latRad) * std::cos(angularDist) +
                                     std::cos(latRad) * std::sin(angularDist) * std::cos(bearingRad));
        
        double outLonRad = lonRad + std::atan2(
            std::sin(bearingRad) * std::sin(angularDist) * std::cos(latRad),
            std::cos(angularDist) - std::sin(latRad) * std::sin(outLatRad));
        
        outLat = radiansToDegrees(outLatRad);
        outLon = radiansToDegrees(outLonRad);
    }


}  // namespace geohelpers

