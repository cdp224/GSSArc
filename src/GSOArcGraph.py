import numpy as np
import matplotlib.pyplot as plt
from Sphere import calc_az_elev_Coordinates

def angle(v1, v2):
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return angle

def calc_gso_position_sph(obs_latitude_deg, sat_longitude_deg):

    # GS Satellite Orbit in km
    Rs = 42164 
    # Earth radius in km
    Re = 6371

    obs_latitude_rad = obs_latitude_deg * np.pi / 180
    sat_longitude_rad = sat_longitude_deg * np.pi / 180

    #Define center of earth
    earth_center = [0,0,0]

    #Define Observer's position
    x_E = 0
    y_E = Re * np.cos(obs_latitude_rad)
    z_E = Re * np.sin(obs_latitude_rad)
    observer = np.array([x_E, y_E, z_E]) #Position of observer
    
    print("Observer latitude: " +  str(obs_latitude_deg) + "°")        
    print("===============================")
    
    print("Satellite longitude: " +  str(np.degrees(sat_longitude_rad)) + "°")        
    #Set position of GSO satellite over a given relative longitude from the observer
    x_S = Rs * np.sin(sat_longitude_rad)
    y_S = Rs * np.cos(sat_longitude_rad)
    z_S = 0
    satellite = np.array([x_S, y_S, z_S])

    print("Satellite overheadposition: " +  str(0) + "°")        
    
    #Set position of satellite over same longitude as observer (i.e., relative longitude 0)
    x_S_eq = Rs * np.sin(0)
    y_S_eq = Rs * np.cos(0)
    z_S_eq = 0
    satellite_eq = np.array([x_S_eq, y_S_eq, z_S_eq])
    
    # The satellite is directly south , so its azimuth is 180 degrees
    ref_long_deg = 180

    # The zenit is when the observer looks up (in opposite direction of earth's center)
    zenit = np.add(observer, np.subtract(observer, earth_center))

    azimuth_deg, elevation_deg  = calc_az_elev_Coordinates(satellite, observer, zenit, satellite_eq, ref_long_deg)
    
    return azimuth_deg, elevation_deg

def calc_gso_position_sg(obs_latitude_deg, sat_longitude_deg):

    # GS Satellite Orbit in km
    Rs = 42164 
    # Earth radius in km
    Re = 6371

    obs_latitude_rad = obs_latitude_deg * np.pi / 180
    sat_longitude_rad = sat_longitude_deg * np.pi / 180

    x_E = 0
    y_E = Re * np.cos(obs_latitude_rad)
    z_E = Re * np.sin(obs_latitude_rad)
    observer = np.array([x_E, y_E, z_E]) #Position of observer
    print("Observer latitude: " +  str(obs_latitude_deg) + "°")        
    print("===============================")
    
    print("Satellite longitude: " +  str(np.degrees(sat_longitude_rad)) + "°")        
    #Set position of satellite
    x_S = Rs * np.sin(sat_longitude_rad)
    y_S = Rs * np.cos(sat_longitude_rad)
    z_S = 0
    satellite = np.array([x_S, y_S, z_S])

    print("Satellite overheadposition: " +  str(0) + "°")        
    #Set position of satellite
    x_S_eq = Rs * np.sin(0)
    y_S_eq = Rs * np.cos(0)
    z_S_eq = 0
    satellite_aq = np.array([x_S_eq, y_S_eq, z_S_eq])

    #put satellites on same sphere around observer than satellite_aq
    vec_observer_satellite = np.subtract(satellite, observer)
    vec_observer_satellite = np.multiply(1/np.linalg.norm(vec_observer_satellite), vec_observer_satellite)

    vec_observer_satellite_aq = np.subtract(satellite_aq, observer)
    vec_observer_satellite_aq = np.multiply(1/np.linalg.norm(vec_observer_satellite_aq), vec_observer_satellite_aq)

    #calculate the zenit: straight over the observer on the sphere.
    vec_observer_zenit = observer
    vec_observer_zenit = np.multiply(1/np.linalg.norm(vec_observer_zenit), vec_observer_zenit)

    elevation_deg = 90-np.degrees(angle(vec_observer_zenit, vec_observer_satellite))
    azimuth_deg=180+np.sign(sat_longitude_rad)*np.degrees(angle(vec_observer_satellite_aq, vec_observer_satellite))

    return azimuth_deg, elevation_deg


def calc_gso_position(obs_latitude_deg, sat_longitude_deg):

    # GS Satellite Orbit in km
    Rs = 42164 
    # Earth radius in km
    Re = 6371

    obs_latitude_rad = obs_latitude_deg * np.pi / 180
    sat_longitude_rad = sat_longitude_deg * np.pi / 180

    x_E = 0
    y_E = Re * np.cos(obs_latitude_rad)
    z_E = Re * np.sin(obs_latitude_rad)
    observer = np.array([x_E, y_E, z_E]) #Position of observer
    print("Observer latitude: " +  str(obs_latitude_deg) + "°")        
    print("===============================")
    
    print("Satellite longitude: " +  str(np.degrees(sat_longitude_rad)) + "°")        
    #Set position of satellite
    x_S = Rs * np.sin(sat_longitude_rad)
    y_S = Rs * np.cos(sat_longitude_rad)
    z_S = 0
    satellite = np.array([x_S, y_S, z_S])

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
    azimuth_deg = np.sign(sat_longitude_rad) * np.degrees(azimuth_rad) + 180
    print("  Satellite azimuth:   " + str(azimuth_deg) + "°")

    return azimuth_deg, elevation_deg

def geo_arc(latitude_deg):
    """
    Simplified computation of azimuth and elevation of possibly visible geostationary satellites
    from a given observer latitude.

    Parameters:
        lat_deg (float): Latitude in degrees
    
    Returns:
        geo_arc (ndarray): Contains azimuth and elevation in degrees
    """

    #Observers latitutde # -90 ... 0 ... 90 degrees (south to north)
    #Observers longitude is 0 #Longitude Range: -180 ... 0 ... 180 degrees (West ... 0 ... East)

    sat_longitude_deg = np.linspace(-90, 90, 100)

    geoarc_a = np.array([[0,0]])

    #iterate over all given satellite longitudes
    for longitude_deg in sat_longitude_deg:

        azimuth_deg, elevation_deg = calc_gso_position_sph(latitude_deg, longitude_deg)
        
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
    plt.rc('font', size=15)
    hfont = {'fontname':'Arial'}
    plt.xlabel('Azimuth [°]', **hfont, size=15)
    plt.ylabel('Elevation [°]', **hfont, size=15)
    plt.title('Visible Geostationary Satellite Arc by Latitude', **hfont, size=20)
    plt.grid(True, linestyle=':')
    L = plt.legend()
    plt.setp(L.texts, family='Arial')
    plt.legend()
    plt.xlim(100, 260)
    plt.xticks(np.arange(105, 256, step=15))
    plt.ylim(0, 45)
    plt.tight_layout()
    print("Save figure to pdf and png into ./output folder.")
    print("The svg-File can be easily dropped into a Word document.")
    plt.savefig('../output/GEOArc.pdf')
    plt.savefig('../output/GEOArc.png')
    plt.savefig('../output/GEOArc.eps')
    plt.savefig('../output/GEOArc.svg')
    print("Open figure window. Close it to terminate.")
    plt.show()


if __name__ == "__main__":
    latitudes = [35, 40, 50, 60, 70, 80]
    plot_geo_arcs(latitudes)
