import Rhino.Geometry as rg
import copy

gap = 5
#currentLevel = 0
#for each crv divide equidistannt
#add a gap in the beginning?? (later)

#distance that triggers the combination of two curves into one
dist_tolerance *= y

crvs = []
crvs.append(crv1)
crvs.append(crv2)

current_in_crv = crv1
current_out_crv = crv2

def CalculateISF(st_lb, st_hb,levels, currentLevel):
    trel = currentLevel/levels
    return st_lb * ( 1/st_lb * st_hb ) ** trel


#function that provides the new points 
# list of points extruded 
def CreatePointsOnDome(prev_points, prev_planes, curve, st_lb, st_hb, levels, currentLevel, width, isInside):
    #we take the points, tangent and planes from prev height
    #Declare ISF for this level:
    iSF = CalculateISF(st_lb, st_hb,levels, currentLevel)
    iSF *= width
    #move based on the tangents
    i = 0
    new_points = []
    adj_new_pts = []
    new_planes = []
    new_params = []
    
    for i in range(len(prev_points)):
        #get the right direction vector
        unitizedDir = prev_planes[i].XAxis * -1
        actualDir = unitizedDir * iSF
        if isInside:
            actualDir *= reductionVal
        new_point = actualDir + prev_points[i]
        new_point.Z += z
        new_points.append(new_point)
        
    #create the Curve
    new_crv = rg.Curve.CreateInterpolatedCurve(new_points, 3)
    #get the tangents and create the new planes
    #redundant but check !!!
    adj_new_pts = new_crv.DivideEquidistant(x/2+gap/2)
    for pt in adj_new_pts:
        output = new_crv.ClosestPoint(pt,0.0001)
        if not None:
            new_params.append(output[1])
        
    for i in new_params:
        p = new_crv.PointAt(i)
        t = new_crv.TangentAt(i)
        v = rg.Vector3d.CrossProduct(t, rg.Vector3d.ZAxis)
        pn =  rg.Plane(p,v,t)
        new_planes.append(pn)
    
    return adj_new_pts, new_params, new_planes, new_crv

#def CheckIntersection(points_in_lvl, points_out_lvl, curve_in, curve_out, st_lb, st_hb, levels, currentLevel, width, dist_tolerance isInside):
#    checked_pts_in = []
#    checked_pts_out = []
#    
#    for pt in points_in_lvl
#        if curve_out.ClosestPoint(pt, dist_tolerance):
#            checked_pts_in.append(pt)

def CheckIntersection(pt, curve_out):
    print curve_out.ClosestPoint(pt, dist_tolerance)[0]
    return curve_out.ClosestPoint(pt, dist_tolerance)[0]


planes_out = []
params_out = []
points_out = []
curves_out = []

planes_in = []
params_in = []
points_in = []
curves_in = []

for crv in crvs:
    params = []
    planes = []
    points = crv.DivideEquidistant(x/2+gap/2)
    for pt in points:
        output = crv.ClosestPoint(pt,0.0001)
        if not None:
            params.append(output[1])
        
    for i in params:
        p = crv.PointAt(i)
        t = crv.TangentAt(i)

        v = rg.Vector3d.CrossProduct(t, rg.Vector3d.ZAxis)
        
        pn =  rg.Plane(p,v,t)

        planes.append(pn)
    
    if crvs.index(crv) == 0:
        planes_in.append(planes)
        params_in.append(params)
        points_in.append(points)
        curves_in.append(crv)
    
    else:
        planes_out.append(planes)
        params_out.append(params)
        points_out.append(points)
        curves_out.append(crv)

#increment the current height
#currentLevel +=1

#take each curve
#iterate of height 



#for outside Curve:
for level in range(levels-1):
    new_pts, new_params, new_planes, new_crv = CreatePointsOnDome(points_out[-1], planes_out[-1], current_out_crv, st_lb, st_hb, levels, level+1,y, False)
    points_out.append(new_pts)
    params_out.append(new_params)
    planes_out.append(new_planes)
    current_out_crv = new_crv
    curves_out.append(new_crv)
    
#for inside Curve:
for level in range(levels-1):
    new_pts, new_params, new_planes, new_crv = CreatePointsOnDome(points_in[-1], planes_in[-1], current_in_crv, st_lb, st_hb, levels, level+1,y ,True)
    points_in.append(new_pts)
    params_in.append(new_params)
    planes_in.append(new_planes)
    current_in_crv = new_crv
    curves_in.append(new_crv)


ready_pts_out = []
ready_params_out = []
ready_planes_out = []
#ready_curves_out = []


#for outer curve 
for i in range(len(points_out)):
    #get a list of point based on i
    pts = [pt for pt in points_out[i][i%2::2]]
    #get a list of params based on i
    params = [param for param in params_out[i][i%2::2]]
    #get a list of planes based on i 
    planes = [plane for plane in planes_out[i][i%2::2]]
    
    ready_pts_out.append(pts)
    ready_params_out.append(params)
    ready_planes_out.append(planes)


ready_pts_in = []
ready_params_in = []
ready_planes_in = []
#ready_curves_in = []

#for inner curve 
for i in range(len(points_in)):
    pts = []
    planes = []
    params = []
    
    #get a list of point based on i
    #print points_in[i]
    points_in[i] = [pt for pt in points_in[i][i%2::2]]
    params_in[i] = [param for param in params_in[i][i%2::2]]
    planes_in[i] = [plane for plane in planes_in[i][i%2::2]]

    for j in range(len(points_in[i])):
        if not CheckIntersection(points_in[i][j], curves_out[i]):
            pts.append(points_in[i][j])
            params.append(params_in[i][j])
            planes.append(planes_in[i][j])
            print "we are here"
            
    ready_pts_in.append(pts)
    ready_params_in.append(params)
    ready_planes_in.append(planes)
