import numpy as np
from scipy import stats
import time

class Euler:

  def __init__(self,n=1,deg=False,**kwargs):

    self.angles=self.__random_euler(n,deg,**kwargs)

  def __random_euler(self,n,deg,**kwargs):
    
    M = stats.special_ortho_group.rvs(dim=3,size=n,**kwargs)
    
    if n==1:
      return np.array([self.__mat2euler(M)])
    else:
      dum=[self.__mat2euler(M[i,:,:]) for i in range(0,n)]
      return np.stack(dum,axis=0)

  def __mat2euler(self,M,deg=False):

    Phi = np.arccos(M[2,2])
    if Phi == 0.0:
      phi1 = np.arctan2(-M[1,0], M[0,0])
      phi2 = 0.0  
    elif Phi == np.pi:
      phi1 = np.arctan2(M[1,0], M[0,0])
      phi2 = 0.0
    else:
      phi1 = np.arctan2(M[2,0], -M[2,1])
      phi2 = np.arctan2(M[0,2], M[1,2])
    
    if phi1 < 0.0:
      phi1 += 2*np.pi
    if phi2 < 0.0:
      phi2 += 2*np.pi
    
    if deg:
      phi1=phi1/np.pi*180.0
      phi2=phi2/np.pi*180.0
      Phi=Phi/np.pi*180

    return np.array([phi1,Phi,phi2])


class TimerError(Exception):
    
    pass
    
class Timer:

    def __init__(self):
        
        self._start_time = None

    def start(self):

        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self,msg='Elapsed time'):

        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        print(f"\n{msg}: {elapsed_time:0.4f} seconds")
    
    def report(self):

        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        return time.perf_counter() - self._start_time
