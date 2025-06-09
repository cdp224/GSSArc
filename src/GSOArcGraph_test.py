# test_my_module.py
from GSOArcGraph import calc_gso_position
from GSOArcGraph import calc_gso_position_sg
from GSOArcGraph import calc_gso_position_sph

def test_both():
    print("Compare both functions")
    astra_long = 19.2
    #From Berne (47N 7.41E) Astra 19.2 longitude
    print("Berne (calc_gso_position)")
    az, el = calc_gso_position(47, astra_long-7.41)
    print("az=" + str(az))
    print("el=" + str(el))
    print("Berne (calc_gso_position_sg)")
    az, el = calc_gso_position_sg(47, astra_long-7.41)
    print("sg_az=" + str(az))
    print("sg_el=" + str(el))
    print("Berne (calc_gso_position_sph)")
    az, el = calc_gso_position_sph(47, astra_long-7.41)
    print("sph_az=" + str(az))
    print("sph_el=" + str(el))

    #From Helsinki (60.1742129N,24.9408216E) Astra 19.2 longitude
    print("Helsinki (calc_gso_position)")
    az, el = calc_gso_position(60.1742129, astra_long-24.9408216)
    print("az=" + str(az))
    print("el=" + str(el))
    print("Helsinki (calc_gso_position_sg)")
    az, el = calc_gso_position_sg(60.1742129, astra_long-24.9408216)
    print("sp_az=" + str(az))
    print("sp_el=" + str(el))
    print("Helsinki (calc_gso_position_sph)")
    az, el = calc_gso_position_sph(60.1742129, astra_long-24.9408216)
    print("sph_az=" + str(az))
    print("sph_el=" + str(el))

    #From Budapest (47.4951479N,19.0572549 E,) Astra 19.2 longitude
    print("Budapest (calc_gso_position)")
    az, el = calc_gso_position(47.4951479, astra_long-19.0572549)
    print("az=" + str(az))
    print("el=" + str(el))
    print("Budapest (calc_gso_position_sg)")
    az, el = calc_gso_position_sg(47.4951479, astra_long-19.0572549)
    print("sg_az=" + str(az))
    print("sg_el=" + str(el))
    print("Budapest (calc_gso_position_sph)")
    az, el = calc_gso_position_sph(47.4951479, astra_long-19.0572549)
    print("sph_az=" + str(az))
    print("sph_el=" + str(el))



def test_calc_gso_position():
    print("test_calc_gso_position")
    astra_long = 19.2
    #From Berne (47N 7.41E) Astra 19.2 longitude
    print("Berne")
    az, el = calc_gso_position(47, astra_long-7.41)
    print("az=" + str(az))
    print("el=" + str(el))
    
    #From Stockholm (59.3179597N,18.0931619E) Astra 19.2 longitude
    print("Stockholm")
    az, el = calc_gso_position(59.3, astra_long-18.09)
    print("az=" + str(az))
    print("el=" + str(el))

    #From Helsinki (60.1742129E,24.9408216E) Astra 19.2 longitude
    print("Helsinki")
    az, el = calc_gso_position(60.1742129, astra_long-24.9408216)
    print("az=" + str(az))
    print("el=" + str(el))

    #assert add(2, 3) == 5
    #assert add(-1, 1) == 0
    #assert add(0, 0) == 0

def test_calc_gso_position_sg():
    print("test_calc_gso_position_sg")
    astra_long = 19.2
    #From Berne (47N 7.41E) Astra 19.2 longitude
    print("Berne")
    az, el = calc_gso_position_sg(47, astra_long-7.41)
    print("az=" + str(az))
    print("el=" + str(el))
    
    #From Stockholm (59.3179597N,18.0931619E) Astra 19.2 longitude
    print("Stockholm")
    az, el = calc_gso_position_sg(59.3, astra_long-18.09)
    print("az=" + str(az))
    print("el=" + str(el))

    #From Helsinki (60.1742129E,24.9408216E) Astra 19.2 longitude
    print("Helsinki")
    az, el = calc_gso_position_sg(60.1742129, astra_long-24.9408216)
    print("az=" + str(az))
    print("el=" + str(el))

    #assert add(2, 3) == 5
    #assert add(-1, 1) == 0
    #assert add(0, 0) == 0

test_calc_gso_position()
test_calc_gso_position_sg()
test_both()
print("All tests passed.")
