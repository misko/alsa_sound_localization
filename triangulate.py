import numpy as np
#import matplotlib.pyplot as plt
from scipy.optimize import fmin,fmin_cg,fmin_bfgs
from scipy.spatial.distance import pdist,squareform



class Triangulator(object):
    def __init__(self,points):
        self.points = points

    def distances(self,point):
        return np.linalg.norm((self.points - point),axis=1)


    def position(self,r,init=None):
        if init is None:
            init = np.average(self.points,axis=0)

        energy = lambda x: np.linalg.norm(self.distances(x)**2-r**2)**2
        grad = lambda x: np.linalg
        x = fmin(energy,init)
        '''
        print x
        print "Dim errors:",(self.distances(x)-r)
        print
        x = fmin_cg(energy,init)
        print x
        print "Dim errors:",(self.distances(x)-r)
        print
        x = fmin_bfgs(energy,init)
        print x
        print "Dim errors:",(self.distances(x)-r)
        print
        '''
        return x


def where_am_i(ps,x1,x2):
    #make an oracle for each point
    os=[]
    ang_dists=np.zeros(ps.shape[0])
    for i in xrange(ps.shape[0]):
        os.append(make_oracle(ps[i,:]))
	ang_dists[i]=os[-1](x1,x2)

    def optim_func(x):
        #x is two points
        ret = 0
	for i in xrange(ps.shape[0]):
            ret+=(ang_dists[i]-os[i](x[0:2],x[2:]))**2
        return ret
   
    init = np.zeros(4)
    x = fmin(optim_func,init)
    print "X1 solved:",x[:2],"actual",x1
    print "X2 solved:",x[2:],"actual",x2
    return x
    

def calibration(oracle):
    points = .2 * np.array([(1,0),(0,1),(0,0),(1,1)]) + .4

    #depend on some magic thing that gives relative distances
    ang_dists = np.zeros((points.shape[0],points.shape[0]))
    for i in xrange(points.shape[0]):
        for j in xrange(i,points.shape[0]):
            ang_dists[i,j] = oracle(points[i,:],points[j,:])
    
    def optim_func(p):
        o = make_oracle(p)
        ret = 0
        for i in xrange(points.shape[0]):
            for j in xrange(i,points.shape[0]):
                ret = ret + (ang_dists[i,j] - o(points[i,:],points[j,:]))**2
        return ret


    init = np.zeros(points.shape[1]+1)
    init[:-1] = np.average(points,axis=0)
    init[-1] = 1
    x = fmin(optim_func,init)
    if x[-1] < 0:
        x[-1] = - x[-1]
    return x
    


def make_oracle(point):
	
    #distance from 2D plane coordinates X to 3D point p    
    def upgrade_norm(x,p):
        return (np.linalg.norm(x-p[:-1])**2 + p[-1]**2)**.5

    # the change in distance to p while moving from X->Y
    return lambda x,y: upgrade_norm(y,point) - upgrade_norm(x,point)
    #return lambda x,y: np.linalg.norm(np.vstack((y,np.array([0]))),point) - np.linalg.norm(np.vstack((x,np.array([0]))),point)


def __main__():
    #calibration test using only relative points using a known movement pattern! ( this is how the robot finds where the cameras are in some absolute coordinate frame)
    #this needs to be done for each camera
    res = calibration(make_oracle(np.array([-5,-5,.1])))
    print "Calibration result:"
    print res

    #triangulation using known points (This is how the cameras find the robot at any given time!
    points = np.array([(1,0,0),(0,1,1),(0,0,2),(1,1,0)])
    t = Triangulator(points)
    init = np.array([1.3,.2,0])
    d = t.distances(init)
    print t.position(d)

    ps=np.random.rand(4,3)*10
    x1=np.random.rand(1,2)*10
    x2=np.random.rand(1,2)*10
    where_am_i(ps,x1,x2)

if __name__ == "__main__":
    __main__()
