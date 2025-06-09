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

   print("Axis Vector Length (must be 1):" + str(ux**2+uy**2+uz**2))

   phi = np.radians(phi_deg)

   rotation_matrix = [[np.cos(phi)+ux**2*(1-np.cos(phi)), ux*uy*(1-np.cos(phi))-uz*np.sin(phi), ux*uz*(1-np.cos(phi))+uy*np.sin(phi)],
                      [uy*ux*(1-np.cos(phi))+uz*np.sin(phi), np.cos(phi)+uy**2*(1-np.cos(phi)),  uy*uz*(1-np.cos(phi))-ux*np.sin(phi)],
                      [uz*ux*(1-np.cos(phi))-uy*np.sin(phi),uz*uy*(1-np.cos(phi))+ux*np.sin(phi), np.cos(phi)+uz**2*(1-np.cos(phi))]]

   rotated_vector = np.matmul(vector, rotation_matrix)

   return rotated_vector

def calc_az_elev_Coordinates(point_of_interest, center_point, zenit_point, reference_point, reference_longitude_deg):
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

    #axis that will be used to rotate POI and reference to the horizon
    vect_rotation_axis_reference = normalize_vector(np.cross(vect_center_zenit, vect_center_reference))
    vect_rotation_axis_poi = normalize_vector(np.cross(vect_center_zenit, vect_center_poi))

    #calculate elevation of POI
    elevation_poi_rad = np.pi/2-angle(vect_center_zenit, vect_center_poi)
    elevation_poi_deg = np.degrees(elevation_poi_rad)
    
    # calcualte the elevation of the reference point.
    elevation_reference_rad = np.pi/2-angle(vect_center_zenit, vect_center_reference)
    elevation_reference_deg = np.degrees(elevation_reference_rad)

    #Goal: bring poi to same down to zero
    rotation_angle_deg = - elevation_poi_deg
    print("POI Rotation angle: "+str(rotation_angle_deg))

    vect_center_poi_rot = rotate(vect_center_poi, vect_rotation_axis_poi, rotation_angle_deg)
    print("POI elevation before rotation " +str(vect_center_poi))
    print(elevation_poi_deg)
    print("POI elevation after rotation (should be abt. zero) " +str(vect_center_poi_rot))
    print(90-np.degrees(angle(vect_center_zenit, vect_center_poi_rot)))

    #Goal: bring reference to same down to zero
    rotation_angle_deg = - elevation_reference_deg
    print("Reference Rotation angle: "+str(rotation_angle_deg))
    vect_center_reference_rot  = rotate(vect_center_reference, vect_rotation_axis_reference, rotation_angle_deg)
    print("Reference elevation before rotation")
    print(elevation_reference_deg)
    print("Reference elevation after rotation (should be abt. zero)")
    print(90-np.degrees(angle(vect_center_zenit, vect_center_reference_rot)))


    #determine sign:
    print("tr(np.dot(normalize_vector(np.cross(vect_center_reference_rot, vect_center_poi_rot)), vect_center_zenit)) ddddd.... " + str(-1*np.dot(normalize_vector(np.cross(vect_center_reference_rot, vect_center_poi_rot)), vect_center_zenit)))
    sign = (-1*np.dot(normalize_vector(np.cross(vect_center_reference_rot, vect_center_poi_rot)), vect_center_zenit))
    if str(sign)=="nan":
        print("nan detected")
        sign = 1
    if sign > 0:
        sign = 1
    else:
        sign = -1
    print("sign: " + str(sign))
    print("plain azi: " + str(np.degrees(angle(vect_center_poi_rot, vect_center_reference_rot))))
    azimuth_deg = reference_longitude_deg + sign*np.degrees(angle(vect_center_poi_rot, vect_center_reference_rot))
    azimuth_deg = np.mod(azimuth_deg, 360)

    print("Elev poi: "+str(elevation_poi_deg))
    print("Elev " + str(elevation_reference_deg))

    print("Azi "+str(azimuth_deg))
    return azimuth_deg, elevation_poi_deg

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

    point_of_interest = [0, 1, 0]
    center_point = [0, 0, 0]
    zenit_point = [0, 0, 1]
    reference_point=[-1, 0, 0]
    reference_longitude_deg = 180
    az, el = calc_az_elev_Coordinates(point_of_interest, center_point, zenit_point, reference_point, reference_longitude_deg)
    print("poi: "+str(point_of_interest))
    print("center: "+str(center_point))
    print("zenit: "+str(zenit_point))
    print("reference: "+str(reference_point))
    print("referenzelongi: "+str(reference_longitude_deg))
    print("AzPOI="+str(az))
    print("ElPOI=",str(el))

    point_of_interest = [0, -1, 0]
    center_point = [0, 0, 0]
    zenit_point = [0, 0, 1]
    reference_point=[-1, 0, 0]
    reference_longitude_deg = 180
    az, el = calc_az_elev_Coordinates(point_of_interest, center_point, zenit_point, reference_point, reference_longitude_deg)
    print("poi: "+str(point_of_interest))
    print("center: "+str(center_point))
    print("zenit: "+str(zenit_point))
    print("reference: "+str(reference_point))
    print("referenzelongi: "+str(reference_longitude_deg))
    print("AzPOI="+str(az))
    print("ElPOI=",str(el))

    point_of_interest = [1, 0, 0]
    center_point = [0, 0, 0]
    zenit_point = [0, 0, 1]
    reference_point=[-1, 0, 0]
    reference_longitude_deg = 180
    az, el = calc_az_elev_Coordinates(point_of_interest, center_point, zenit_point, reference_point, reference_longitude_deg)
    print("poi: "+str(point_of_interest))
    print("center: "+str(center_point))
    print("zenit: "+str(zenit_point))
    print("reference: "+str(reference_point))
    print("referenzelongi: "+str(reference_longitude_deg))
    print("AzPOI="+str(az))
    print("ElPOI=",str(el))

    point_of_interest = [1, 0, 1]
    center_point = [0, 0, 0]
    zenit_point = [0, 0, 1]
    reference_point=[-1, 0, 0]
    reference_longitude_deg = 180
    az, el = calc_az_elev_Coordinates(point_of_interest, center_point, zenit_point, reference_point, reference_longitude_deg)
    print("poi: "+str(point_of_interest))
    print("center: "+str(center_point))
    print("zenit: "+str(zenit_point))
    print("reference: "+str(reference_point))
    print("referenzelongi: "+str(reference_longitude_deg))
    print("AzPOI="+str(az))
    print("ElPOI=",str(el))


    point_of_interest = [1, 0, -200]
    center_point = [0, 0, 0]
    zenit_point = [0, 0, 1]
    reference_point=[-1, 0, 0]
    reference_longitude_deg = 180
    az, el = calc_az_elev_Coordinates(point_of_interest, center_point, zenit_point, reference_point, reference_longitude_deg)
    print("poi: "+str(point_of_interest))
    print("center: "+str(center_point))
    print("zenit: "+str(zenit_point))
    print("reference: "+str(reference_point))
    print("referenzelongi: "+str(reference_longitude_deg))
    print("AzPOI="+str(az))
    print("ElPOI=",str(el))


    