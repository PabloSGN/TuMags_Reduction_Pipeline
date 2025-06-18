from scipy.interpolate import interp1d
import numpy as np

def step(original_curve,center = None,amplitude = None):
        # Up step: amplitude is percent, applied to all x >= center
    modified_curve=np.copy(original_curve)
    for i in range(len(center)):
        modified_curve[int(center[i]):] = modified_curve[int(center[i]):] + amplitude[i] 
    return modified_curve


def modified_curve(center = [1539,1540],amplitude = [1.0,1.0]):
    MAX_VALUE = 4095
    x = np.arange(0, MAX_VALUE + 1)
    modified_curve = step(x,center=center,amplitude=amplitude)
    # Create the interpolation function
    modify_linear = interp1d(x, modified_curve, kind='cubic', bounds_error=False, fill_value=(modified_curve[0], modified_curve[-1]))
    return modify_linear

def modify_curve(modify_linear,input_values):
    return modify_linear(input_values)