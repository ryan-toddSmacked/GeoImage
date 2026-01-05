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
    // GeoTIFF Tag IDs
    //--------------------------------------------------------------------------

    constexpr uint16_t TIFFTAG_GEOKEYDIRECTORY = 34735;   // GeoKeyDirectoryTag
    constexpr uint16_t TIFFTAG_GEODOUBLEPARAMS = 34736;   // GeoDoubleParamsTag
    constexpr uint16_t TIFFTAG_GEOASCIIPARAMS = 34737;    // GeoAsciiParamsTag
    constexpr uint16_t TIFFTAG_MODELPIXELSCALE = 33550;   // ModelPixelScaleTag
    constexpr uint16_t TIFFTAG_MODELTIEPOINT = 33922;     // ModelTiepointTag
    constexpr uint16_t TIFFTAG_MODELTRANSFORMATION = 34264; // ModelTransformationTag

    //--------------------------------------------------------------------------
    // GeoKey IDs
    //--------------------------------------------------------------------------

    constexpr uint16_t GTModelTypeGeoKey = 1024;
    constexpr uint16_t GTRasterTypeGeoKey = 1025;
    constexpr uint16_t GTCitationGeoKey = 1026;
    constexpr uint16_t GeographicTypeGeoKey = 2048;
    constexpr uint16_t GeogCitationGeoKey = 2049;
    constexpr uint16_t GeogGeodeticDatumGeoKey = 2050;
    constexpr uint16_t GeogPrimeMeridianGeoKey = 2051;
    constexpr uint16_t GeogAngularUnitsGeoKey = 2054;
    constexpr uint16_t GeogEllipsoidGeoKey = 2056;
    constexpr uint16_t GeogSemiMajorAxisGeoKey = 2057;
    constexpr uint16_t GeogSemiMinorAxisGeoKey = 2058;
    constexpr uint16_t GeogInvFlatteningGeoKey = 2059;
    constexpr uint16_t ProjectedCSTypeGeoKey = 3072;
    constexpr uint16_t ProjLinearUnitsGeoKey = 3076;

    //--------------------------------------------------------------------------
    // GeoKey Values
    //--------------------------------------------------------------------------

    // Model types (GTModelTypeGeoKey)
    constexpr uint16_t ModelTypeProjected = 1;
    constexpr uint16_t ModelTypeGeographic = 2;
    constexpr uint16_t ModelTypeGeocentric = 3;

    // Raster types (GTRasterTypeGeoKey)
    constexpr uint16_t RasterPixelIsArea = 1;
    constexpr uint16_t RasterPixelIsPoint = 2;

    // Angular units (GeogAngularUnitsGeoKey)
    constexpr uint16_t Angular_Radian = 9101;
    constexpr uint16_t Angular_Degree = 9102;
    constexpr uint16_t Angular_Arc_Minute = 9103;
    constexpr uint16_t Angular_Arc_Second = 9104;

    // Linear units (ProjLinearUnitsGeoKey)
    constexpr uint16_t Linear_Meter = 9001;
    constexpr uint16_t Linear_Foot = 9002;
    constexpr uint16_t Linear_Foot_US_Survey = 9003;

    //--------------------------------------------------------------------------
    // EPSG:4326 (WGS84 Geographic) Constants
    //--------------------------------------------------------------------------

    constexpr uint16_t EPSG_4326 = 4326;                           // Geographic CRS code
    constexpr uint16_t EPSG_4326_DATUM = 6326;                     // WGS84 datum code
    constexpr uint16_t EPSG_4326_ELLIPSOID = 7030;                 // WGS84 ellipsoid code
    constexpr uint16_t EPSG_4326_PRIME_MERIDIAN = 8901;            // Greenwich prime meridian

    // GeoDoubleParams values for EPSG:4326 (WGS84)
    constexpr double EPSG_4326_SEMI_MAJOR_AXIS = 6378137.0;        // Semi-major axis in meters
    constexpr double EPSG_4326_SEMI_MINOR_AXIS = 6356752.314245;   // Semi-minor axis in meters
    constexpr double EPSG_4326_INV_FLATTENING = 298.257223563;     // Inverse flattening

    // Pre-built GeoKeyDirectory for EPSG:4326 (as array of uint16_t values)
    // Format: [KeyDirectoryVersion, KeyRevision, MinorRevision, NumberOfKeys, 
    //          Key1, TIFFTagLocation1, Count1, Value1, ...]
    constexpr std::array<uint16_t, 24> EPSG_4326_GeoKeyDirectory = {
        1, 1, 0, 5,                                    // Header: version 1.1.0, 5 keys
        GTModelTypeGeoKey, 0, 1, ModelTypeGeographic,  // Model = Geographic
        GTRasterTypeGeoKey, 0, 1, RasterPixelIsArea,   // Raster = PixelIsArea
        GeographicTypeGeoKey, 0, 1, EPSG_4326,         // Geographic CRS = EPSG:4326
        GeogGeodeticDatumGeoKey, 0, 1, EPSG_4326_DATUM,// Datum = WGS84
        GeogAngularUnitsGeoKey, 0, 1, Angular_Degree   // Units = Degrees
    };

    // Pre-built GeoDoubleParams for EPSG:4326 (if explicit ellipsoid params needed)
    constexpr std::array<double, 3> EPSG_4326_GeoDoubleParams = {
        EPSG_4326_SEMI_MAJOR_AXIS,
        EPSG_4326_SEMI_MINOR_AXIS,
        EPSG_4326_INV_FLATTENING
    };

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

