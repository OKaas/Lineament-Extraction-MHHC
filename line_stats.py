# =====================================================================================
# | Original MHHC algorithm for spatial clustering of the line segments
# | as a tool for lineament extraction.
# |
# | Jakub Silhavy, 2016
# | Faculty of Applied Sciences, University of West Bohemia, Pilsen, Czech Republic
# =====================================================================================

# function to compute statistic of two line segments
import math
class point:
  def __init__(self, X, Y):
    self.X = X
    self.Y = Y
    
def lineLength(A,B):
  return math.sqrt(math.pow(A.X-B.X,2)+math.pow(A.Y-B.Y,2))
  
def lineAzimuth(A,B):
  dY = (A.Y - B.Y)
  dX = (A.X - B.X)
  atan2 = math.atan2(dY,dX)
  if (atan2<0):
    atan2+=math.pi
  if (atan2<=math.pi/2):
    atan2=math.pi/2-atan2
  if (atan2>math.pi/2):
    atan2=3*math.pi/2-atan2
  return math.degrees(atan2)

# orthogonal distance between point C and line segment A,B
def orthoDist(C, A, B):
  a = lineLength(B,C)
  b = lineLength(A,C)
  c = lineLength(A,B)
  s = (a+b+c)/2
  sqrtDomain = s*(s-a)*(s-b)*(s-c)
  if (sqrtDomain > 0):
    S = math.sqrt(s*(s-a)*(s-b)*(s-c))
  else:
    S = 0
  return 2*S/c

def lineStat(A,B,C,D):
  # postupne kontrola parametru - jeden vyjde nula - ostatni nepocitam!
  # lenght ratio
  lenghtAB = lineLength(A,B)
  lenghtCD = lineLength(C,D)
  
  # ratio 
  dL = (lenghtAB / lenghtCD)
  if dL > 1:
    dL = 1/dL

  # azimuth difference
  sigmaMin = 30*math.pi/180
  sigmaAB = lineAzimuth(A,B)
  sigmaCD = lineAzimuth(C,D)
  dSigma = math.fabs(sigmaAB - sigmaCD)
  sigmaRatio =  1 - dSigma/sigmaMin

  # orthogonal distance
  orthoA = orthoDist(A, C, D)
  orthoB = orthoDist(B, C, D)
  orthoC = orthoDist(C, A, B)
  orthoD = orthoDist(D, A, B)
  
  orthoAv = (orthoA+orthoB+orthoC+orthoD)/4
  scaleNumber = 50 # 1:50 000
  mmTolerance = 4 # mm
  orthoRatio = 1 - orthoAv/(scaleNumber*mmTolerance)

  # vahy
  angleW = 0.5
  distW = 0.25
  lengthW = 0.25
  
  return (angleW*sigmaRatio+distW*orthoRatio+lengthW*dL)/(angleW+distW+lengthW)