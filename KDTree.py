import math,random

# return list of (item,cnt) sorted by counts

def popularity(lst):
  hash = {}
  for x in lst:
    if x not in hash: hash[x] = 0
    hash[x] += 1
  data = [(hash[x],x) for x in hash.keys()]
  data.sort(reverse=True)
  data = [(y,x) for (x,y) in data]
  return data

def frand_pm1(): return 2.*random.random()-1.

def dist(f,s):
  a = f[0] 
  b = f[1] 
  c = f[2] 
  d = s[0] 
  e = s[1] 
  f = s[2] 
# def dist([a,b,c],[d,e,f]): 
  return math.sqrt(distsq([a,b,c],[d,e,f]))

def distsq(f,s):
  a = f[0] 
  b = f[1] 
  c = f[2] 
  d = s[0] 
  e = s[1] 
  f = s[2] 
# def distsq((a,b,c),(d,e,f)):
  return (a-d)*(a-d)+(b-e)*(b-e)+(c-f)*(c-f)

def remove_duplicates(list):
  unique = []
  for x in list:
    if not x in unique: unique.append(x)
  return unique

def intersection(a,b):
  subset = []
  for i in a:
    if i in b: subset.append(i)
  return subset

def pdb_rec(anum,anam,rnam,rnum,chn,x,y,z,occ,bf):
  return "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" % \
      (anum,anam,rnam,chn,rnum%10000,x,y,z,occ,bf)
    
def diff(a,b): return [a[k]-b[k] for k in xrange(3)]

def vsum(a,b): return [a[k]+b[k] for k in xrange(3)]
  
def prod(mat,v):
  w = [0 for i in xrange(3)]
  for i in xrange(3):
    for j in xrange(3): w[i] += mat[i][j]*v[j]
  return w

def cross(u,v):
  a = u[1]*v[2]-u[2]*v[1]
  b = u[2]*v[0]-u[0]*v[2]
  c = u[0]*v[1]-u[1]*v[0]
  return [a,b,c]

def my_dot(a,b): return sum([x*y for x,y in zip(a,b)])

def mult_mat_vec(M,v):
  w = []
  for i in range(3):
    z = 0
    for j in range(3): z += M[i][j]*v[j]
    w.append(z)
  return w

def angle(u,v):
  a2 = my_dot(u,u)
  b2 = my_dot(v,v)
  c2 = my_dot(u,v)
  if a2<=0 or b2<=0: return 0
  d = c2/math.sqrt(a2*b2)
  if d<-1: d = -1
  if d>1: d = 1
  d = math.acos(d)
  return 180*d/3.14159

def angle3(a,b,c):
  u = diff(a,b)
  v = diff(c,b)
  return angle(u,v)

# from p. 75 in Altmann

def rotation_matrix(dir,ang): # in degrees

  (nx,ny,nz) = (dir[0],dir[1],dir[2])
  mag = math.sqrt(nx*nx+ny*ny+nz*nz)
  (nx,ny,nz) = (nx/mag,ny/mag,nz/mag) # normalize to unit vector
  ang *= 3.1415927/180.0
  sn = math.sin(ang)
  sq = math.sin(0.5*ang)
  sq *= sq

  m = [[0 for i in xrange(3)] for j in xrange(3)]

  m[0][0] = 1-2*(ny*ny+nz*nz)*sq
  m[0][1] = -nz*sn+2*nx*ny*sq
  m[0][2] = ny*sn+2*nx*nz*sq

  m[1][0] = nz*sn+2*nx*ny*sq
  m[1][1] = 1-2*(nx*nx+nz*nz)*sq
  m[1][2] = -nx*sn+2*ny*nz*sq

  m[2][0] = -ny*sn+2*nx*nz*sq
  m[2][1] = nx*sn+2*ny*nz*sq
  m[2][2] = 1-2*(nx*nx+ny*ny)*sq

  return m

# based on Andrew Moore's tutorial on kd-trees:
#   http://www.autonlab.org/autonweb/14665/version/2/part/5/data/moore-tutorial.pdf?branch=main&language=en

# example of how to use this:
#  kdtree = KDTree([x.co for x in receptor])
#  (dsq,coord) = kdtree.nearest(prot[j].co)
#  j = kdtree.index_of(coord)

class KDTree:

  # members: leaf (0/1), dim, val (float), cutpt (coord)
  #          left, right (KDTrees)
  #          mins,maxs (coords)

  def __init__(self,coords):
    self.INFY = 999999
    self.inverse = {}
    for i in range(len(coords)):
      self.inverse[coords[i][0],coords[i][1],coords[i][2]] = i
    if len(coords)==0: self.leaf = 1; return
    else: self.leaf = 0
    (self.mins,self.maxs) = self.compute_extremes(coords)
    diffs = [x-y for x,y in zip(self.maxs,self.mins)]
    if diffs[0]>=diffs[1] and diffs[0]>=diffs[2]: self.dim = 0
    elif diffs[1]>=diffs[0] and diffs[1]>=diffs[2]: self.dim = 1
    else: self.dim = 2
    self.val = (self.maxs[self.dim]+self.mins[self.dim])/2.0

    best = None
    for co in coords:
      Del = math.fabs(co[self.dim]-self.val)
      if best==None or Del<delbest: (best,delbest) = (co,Del)
    self.cutpt = best
    self.val = self.cutpt[self.dim]
    (lset,rset)=([],[])
    for co in coords:
      if co==best: continue # exclude cutpt from either sub-tree
      if co[self.dim]<=self.val: lset.append(co)
      else: rset.append(co)
    self.left = KDTree(lset)
    self.right = KDTree(rset)

  def num_nodes(self):
    if self.leaf==1: return 1
    return 1+self.left.num_nodes()+self.right.num_nodes()

  def num_leaves(self):
    if self.leaf==1: return 1
    return self.left.num_leaves()+self.right.num_leaves()

  def index_of(self,co): return self.inverse[co[0],co[1],co[2]]

  # returns INFY if no pt found within dmax, else returns distance-squared

  def nearest(self,targ):
    return self.nearest_(targ,self.mins,self.maxs,self.INFY)

  def nearest_(self,targ,mins,maxs,dmax):
    if self.leaf==1: return (self.INFY,None)

    (minN,minF) = (mins[0:3],mins[0:3])
    (maxN,maxF) = (maxs[0:3],maxs[0:3])
    if targ[self.dim]<=self.val: 
      (nearer,further) = (self.left,self.right)
      (maxN[self.dim],minF[self.dim]) = (self.val,self.val)
    else: 
      (nearer,further) = (self.right,self.left)
      (minN[self.dim],maxF[self.dim]) = (self.val,self.val)

    (dsq,output) = nearer.nearest_(targ,minN,maxN,dmax)
    if dsq<dmax: dmax = dsq
    if dmax<(targ[self.dim]-self.val)*(targ[self.dim]-self.val): 
      return (dsq,output)
    d2 = distsq(targ,self.cutpt)
    if d2<dsq:
      output = self.cutpt
      dmax = dsq = d2

    (d3,out) = further.nearest_(targ,minF,maxF,dmax)
    if d3<dsq: (dsq,output) = (d3,out)
    return (dsq,output)

  def near(self,coord):
    (d,pt) = self.nearest(coord)
    if d<DMESH*DMESH: return 1
    return 0

  def compute_extremes(self,coords):
    xs = map(lambda x: x[0],coords)
    ys = map(lambda x: x[1],coords)
    zs = map(lambda x: x[2],coords)
    # print(f"xs is {list(xs)}")
    # print(f"ys is {list(ys)}")
    # print(f"zs is {list(zs)}")

    xlist = list(xs)
    ylist = list(ys)
    zlist = list(zs)

    # print(f"xs is {xlist}")
    # print(f"ys is {(ylist)}")
    # print(f"zs is {(zlist)}")

    mins = [min(xlist),min(ylist),min(zlist)]
    maxs = [max(xlist),max(ylist),max(zlist)]
    return (mins,maxs)

  def range_query_radius(self,co,rad):
    lows = map(lambda x: x-rad,co)
    highs = map(lambda x: x+rad,co)
    hits = self.range_query(lows,highs)
    return filter(lambda x: distsq(x,co)<rad*rad,hits)

  def range_query(self,lows,highs):
    if self.leaf==1: return []
    for i in range(3):
      if lows[i]>self.maxs[i]: return []
      if highs[i]<self.mins[i]: return []
    # there must be some overlap volume
    a = self.left.range_query(lows,highs)
    b = self.right.range_query(lows,highs)
    c = [self.cutpt]
    for i in range(3):
      if self.cutpt[i]<lows[i] or self.cutpt[i]>highs[i]: c = []
    return a+b+c

