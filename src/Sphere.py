import numpy as np

def normalize_vector(vector):
    return np.multiply(1/np.linalg.norm(vector), vector)
    
def angle(v1, v2):
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return angle

def rotate(vector, axis, phi_deg):
   
   ux=axis[0]
   uy=axis[1]
   uz=axis[2]

   phi = np.radians(phi_deg)

   rotation_matrix = [[np.cos(phi)+ux**2*(1-np.cos(phi)), ux*uy*(1-np.cos(phi))-uz*np.sin(phi), ux*uz*(1-np.cos(phi))+uy*np.sin(phi)],
                      [uy*ux*(1-np.cos(phi))+uz*np.sin(phi), np.cos(phi)+uy**2*(1-np.cos(phi)),  uy*uz*(1-np.cos(phi))-ux*np.sin(phi)],
                      [uz*ux*(1-np.cos(phi))-uy*np.sin(phi),uz*uy*(1-np.cos(phi))+ux*np.sin(phi), np.cos(phi)+uz**2*(1-np.cos(phi))]]

   rotated_vector = np.matmul(vector, rotation_matrix)

   return rotated_vector

def az_elev_Coordinates(point_of_interest, center_point, zenit_point, reference_point, reference_longitude_deg):
    # point_of_interest: satellite somewhere
    # center_point: center_point of sphere
    # zenit_point: zenit of sphere
    # reference_point: point for e.g., the azimuth reference (point that points into reference_longitude direction)

    #Vectors:
    vect_center_poi = np.subtract(point_of_interest, center_point)
    vect_center_zenit = np.subtract(zenit_point, center_point)
    vect_center_reference = np.subtract(reference_point, center_point)
    
    vect_center_poi = normalize_vector(vect_center_poi)
    vect_center_zenit = normalize_vector(vect_center_zenit)
    vect_center_reference = normalize_vector(vect_center_reference)
    print("POI")
    print(vect_center_poi)
    print("Zenit")
    print(vect_center_zenit)
    print("Reference")
    print(vect_center_reference)

    vect_west = normalize_vector(np.cross(vect_center_zenit, vect_center_reference))


    elevation_deg = 90-np.degrees(angle(vect_center_zenit, vect_center_poi))
    print("ANgle "+str(np.degrees(angle(vect_center_zenit, vect_center_poi))))
    elevation_reference_deg = 90-np.degrees(angle(vect_center_zenit, vect_center_reference))

    #Goal: bring poi to same elevation in order to calculate the azimuth
    rotation_angle_rad = elevation_reference_deg-elevation_deg
    print("Rotation angle: "+str(rotation_angle_rad))

    vect_center_poi_rot = rotate(vect_center_poi, vect_west, rotation_angle_rad)
    print("elevation before rotation")
    print(elevation_deg)
    
    print("elevation after rotation")
    print(90-np.degrees(angle(vect_center_zenit, vect_center_poi_rot)))

    print("elevation of reference")
    print(elevation_reference_deg)


    azimuth_deg = reference_longitude_deg + np.degrees(angle(vect_center_poi_rot, vect_center_reference))

    print("Elev poi: "+str(elevation_deg))
    print("Elev " + str(elevation_reference_deg))

    print("Azi "+str(azimuth_deg))
    return elevation_deg, azimuth_deg

if __name__ == "__main__":

    vect=[1,0,0]
    axis=[0,0,1]
    phi=90
    print("vec")
    print(vect)
    print("axis")
    print(axis)
    print(phi)
    vect=rotate(vect, axis, phi)
    print(vect)
    print(phi)
    vect=rotate(vect, axis, phi)
    print(vect)
    print(phi)
    vect=rotate(vect, axis, phi)
    print(vect)

    