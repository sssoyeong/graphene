from predict2 import *
import matplotlib.pyplot as plt
import numpy as np

input_list = '61,37,13,21,17,41,45,65,69'
#input_list = '61,60'
atom_category = 'B'
broad = 0.05
erange = 10.
r_cut = 15.

predicted_data, expected_list, prdf_inp, prdf_org = predict(atom_category, input_list, broad, erange)
x_list = np.linspace(-erange,erange,1001)

#print predicted_data
#print len(predicted_data)

fig = plt.figure(1, figsize=(5,7))

fig1 = fig.add_subplot(211)
fig1.plot(x_list, predicted_data, 'k-', label='prediction')
fig1.plot(x_list, expected_list, 'r-', label='expectation')
fig1.legend()
fig1.axis([-10,10,0,20])

fig2 = fig.add_subplot(212)
fig2.plot(np.linspace(0,r_cut,len(prdf_inp)), prdf_inp, 'k-', label='PRDF_input')
fig2.plot(np.linspace(0,r_cut,len(prdf_inp)), prdf_org, 'r-', label='PRDF_orig')
fig2.legend()
#fig2.axis([0,10,0,10])
#fig.savefig('result.png')
plt.show()

