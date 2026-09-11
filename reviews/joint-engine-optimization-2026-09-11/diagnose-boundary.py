import sys
sys.path.insert(0,'/work/tools')
import numpy as np
from benchmark_ten_shifts_missing import fixture
import nwkit.shift_joint_fit as f
from nwkit.shift_native_model import ShiftLayout
from nwkit.shift_native_fit import NativeFitOptions
from scipy.optimize._numdiff import approx_derivative
orig=f.minimize

def observed(fun,x,**kwargs):
 r=orig(fun,x,**kwargs)
 print('OPT',r.success,r.fun,r.x,r.jac,flush=True)
 if kwargs.get('jac') is True:
  for h in [1e-3,1e-4,1e-5,1e-6]:
   lo,hi=zip(*kwargs['bounds']);lo=[-np.inf if a is None else a for a in lo];hi=[np.inf if a is None else a for a in hi]
   g=approx_derivative(lambda v:fun(v)[0],r.x,abs_step=h,bounds=(lo,hi),method='3-point')
   print('GRAD',h,g,flush=True)
 return r
f.minimize=observed
data,_=fixture(dict(truth='shared',replicate=0,traits=2,missing_rate=.2))
layout=ShiftLayout.build(data.tree,[16,17,19,20,23,25,27,29,107,161,178])
print(f._general_fit(data,layout,NativeFitOptions(trait_covariance='full',alpha_model='shared',estimate_measurement_error=True),np.zeros(2),None,np.zeros(2),f.JointFitContext(data,'OUfixedRoot')))
