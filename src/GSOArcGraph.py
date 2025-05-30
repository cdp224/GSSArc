import numpy as np
import matplotlib.pyplot as plt


def angle(v1, v2):
    # v1 is your first vector
    # v2 is your second vector
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return angle

def geo_arc(latitude_deg):
    """
    Simplified computation of azimuth and elevation of possibly visible geostationary satellites
    from a given observer latitude.

    Parameters:
        lat_deg (float): Latitude in degrees
    
    Returns:
        geo_arc (ndarray): Contains azimuth and elevation in degrees
    """
    # GS Satellite Orbit in km
    Rs = 42164 
    # Earth radius in km
    Re = 6371

    #Observers latitutde # -90 ... 0 ... 90 degrees (south to north)
    #Observers longitude is 0 #Longitude Range: -180 ... 0 ... 180 degrees (West ... 0 ... East)

    centerPoint = np.array([0, 0, 0]) # Center point of earth
    #latitude_deg = 59
    obs_latitude_rad = latitude_deg * np.pi/180   
    sat_longitude_rad = np.radians(np.linspace(-90, 90, 100))

    x_E = 0
    y_E = Re * np.cos(obs_latitude_rad)
    z_E = Re * np.sin(obs_latitude_rad)
    observer = np.array([x_E, y_E, z_E]) #Position of observer
    print("Observer latitude: " +  str(latitude_deg) + "°")        
    print("===============================")

    satellite_a = np.array([[0,0,0]])
    geoarc_a = np.array([[0,0]])

    #iterate over all given satellite longitudes
    for longitude_rad in sat_longitude_rad:

        print("Satellite longitude: " +  str(np.degrees(longitude_rad)) + "°")        
        #Set position of satellite
        x_S = Rs * np.sin(longitude_rad)
        y_S = Rs * np.cos(longitude_rad)
        z_S = 0

        satellite = np.array([x_S, y_S, z_S])
        satellite_a = np.concatenate((satellite_a, [satellite]), axis=0)

        #Calculate point on observers horizontal plane where the line center to satellite penetrates 
        r=(y_E*y_E+z_E*z_E)/(y_E*y_S)
        M = np.array(np.multiply(r,[x_S,y_S,0]))

        #calculate elevation (how high is the satellite over the horizon)
        #Angle between two vectors: OS and OM
        vec_observer_satellite = np.subtract(satellite, observer)
        vec_observer_M = np.subtract(M, observer)
        elevation_rad = angle(vec_observer_M, vec_observer_satellite)
        elevation_deg = np.degrees(elevation_rad)
        elevation_sig = np.sign(vec_observer_satellite[1]-vec_observer_M[1])
        elevation_deg = elevation_sig * elevation_deg
        print("  Satellite elevation: " + str(elevation_deg) + "°")

        #calculate the azimuth 
        #Angle between two vectors OMa and OM sat@0 (0 Rs 0)
        M_aequator = np.array([0,(y_E*y_E+z_E*z_E)/(y_E),0])
        vec_observer_M_a = np.subtract(M_aequator, observer)
        azimuth_rad = angle(vec_observer_M_a, vec_observer_M)
        # azimuth + 180 to accomodate for standard direction to south
        azimuth_deg = np.sign(longitude_rad) * np.degrees(azimuth_rad) + 180
        print("  Satellite azimuth:   " + str(azimuth_deg) + "°")
        geoarc_point=np.array([azimuth_deg, elevation_deg])
        geoarc_a = np.concatenate((geoarc_a, [geoarc_point]), axis=0)

    geoarc_a = geoarc_a[1:]
    return geoarc_a
#    return np.degrees(azimuth), np.degrees(elevation)

def plot_geo_arcs(latitudes):

    plt.figure(figsize=(10, 6))

    # Safe line styles for B/W printing
    line_styles = ['-', '--', '-.', ':', '-', '--',]

    for idx, lat in enumerate(latitudes):
        geoarc = geo_arc(lat)
        az=geoarc[:,0]
        el=geoarc[:,1]

        # Sanity check: ensure both are 1D arrays
        if az.ndim != 1 or el.ndim != 1:
            print(f"Invalid output shape for latitude {lat}: az.shape={az.shape}, el.shape={el.shape}")
            continue

        # Filter only valid, positive elevation angles
        valid_mask = np.isfinite(az) & np.isfinite(el) & (el > 0)
        az_valid = np.array(az[valid_mask])
        el_valid = np.array(el[valid_mask])
        
        if az_valid.size == 0 or el_valid.size == 0:
            print(f"No valid data to plot for latitude {lat}°.")
            continue
        
        plt.plot(
            az_valid,
            el_valid,
            line_styles[idx % len(line_styles)],
            label=f'{lat}° N',
            linewidth=1
        )

    plt.xlabel('Azimuth [°]')
    plt.ylabel('Elevation [°]')
    plt.title('Visible Geostationary Satellite Arc by Latitude')
    plt.grid(True, linestyle=':')
    plt.legend()
    plt.xlim(100, 260)
    plt.ylim(0, 45)
    plt.tight_layout()
    print("Save figure to pdf and png into ./output folder.")
    plt.savefig('../output/GEOArc.pdf')
    plt.savefig('../output/GEOArc.png')
    print("Open figure window. Close it to terminate.")
    plt.show()


if __name__ == "__main__":
    latitudes = [35, 40, 50, 60, 70, 80]
    plot_geo_arcs(latitudes)
