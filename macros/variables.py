""" File to especify the variables that are going to be plotted 

>> The structure should appear as follows:
We have a variable dict: variable = {}
variable['variable1_key'] = {'name' : 'name in the file or expression',
                             'range' : [0, 1],
                             'label' : 'label'}



>> NOTE: Expression syntax
Before each variable it is required to write event. such as
event.PV_vx - event.RefittedPV_vx


"""

variables['PV_vx_diff'] = {'name' : 'event.PV_vx - event.RefittedPV_vx',
                      'range' : [-0.02, 0.02],
                      'binning' : 100,
                      'label' : 'Original PV_{x} - Refitted PV_{x}'}


variables['PV_vy_diff'] = {'name' : 'event.PV_vy - event.RefittedPV_vy',
                      'range' : [-0.02, 0.02],
                      'binning' : 100,
                      'label' : 'Original PV_{y} - Refitted PV_{y}'}


variables['PV_vz_diff'] = {'name' : 'event.PV_vz - event.RefittedPV_vz',
                      'range' : [-0.02, 0.02],
                      'binning' : 100,
                      'label' : 'Original PV_{z} - Refitted PV_{z}'}


variables['PV_vT_diff'] = {'name' : 'sqrt(event.PV_vx*event.PV_vx + event.PV_vy*event.PV_vy) - sqrt(event.RefittedPV_vx*event.RefittedPV_vx + event.RefittedPV_vy*event.RefittedPV_vy)',
                      'range' : [-0.05, 0.05],
                      'binning' : 100,
                      'label' : 'Original |PV_{T}| - Refitted |PV_{T}|'}



