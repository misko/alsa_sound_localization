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


def where_am_i(ps,xs):
    #make an oracle for each point
    os=[]
    ang_dists=np.zeros((ps.shape[0],xs.shape[0]-1))
    for i in xrange(ps.shape[0]):
        os.append(make_oracle(ps[i,:]))
        for j in xrange(xs.shape[0]-1):
            ang_dists[i,j]=os[-1](xs[j,:],xs[j+1,:])

    def optim_func(x):
        #x is two points
        ret = 0
        for j in xrange(xs.shape[0]-1):
	    for i in xrange(ps.shape[0]):
                fr=x[j*2:(j+1)*2]
                to=x[(j+1)*2:(j+2)*2]
                ret+=(ang_dists[i,j]-os[i](fr,to))**2
        return ret
   
    init = xs*0
    #init[:2]=x1+5
    #init[2:]=x2
    x = fmin(optim_func,init)
    print x.reshape(xs.shape)
    print xs
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

    ps=(np.random.rand(4,3)-0.5)*10
    xs=(np.random.rand(6,2)-0.5)
    where_am_i(ps,xs)

if __name__ == "__main__":
    __main__()
