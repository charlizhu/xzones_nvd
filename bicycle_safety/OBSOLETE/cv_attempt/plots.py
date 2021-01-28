import statistics
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

toocheap = [12.5]*36+[37.5]*4+[62.5]*2
goodbargain = [12.5]*5+[37.5]*15+[62.5]*17+[87.5]*5
bithigh = [12.5]+[37.5]*6+[62.5]*11+[87.5]*14+[112.5]*4+[137.5]*4+[162.5]+[187.5]
toohigh = [12.5]*2+[62.5]*4+[87.4]*8+[112.5]*3+[137.5]*7+[162.5]*3+[187.5]*4+[212.5]*5+[237.5]+[250]*5

# print(toocheap)
# print(goodbargain)
# print(bithigh)
# print(toohigh)

toocheap_mean = statistics.mean(toocheap)
goodbargain_mean = statistics.mean(goodbargain)
bithigh_mean = statistics.mean(bithigh)
toohigh_mean = statistics.mean(toohigh)

# print(toocheap_mean,goodbargain_mean,bithigh_mean,toohigh_mean)

toocheap_stdev = statistics.stdev(toocheap)
goodbargain_stdev = statistics.stdev(goodbargain)
bithigh_stdev = statistics.stdev(bithigh)
toohigh_stdev = statistics.stdev(toohigh)

# print(toocheap_stdev,goodbargain_stdev,bithigh_stdev,toohigh_stdev)

toocheap_x = np.linspace(toocheap_mean - 3*toocheap_stdev, toocheap_mean + 3*toocheap_stdev, 100)
goodbargain_x = np.linspace(goodbargain_mean-3*goodbargain_stdev, goodbargain_mean+3*goodbargain_stdev, 100)
bithigh_x = np.linspace(bithigh_mean-3*bithigh_stdev,bithigh_mean+3*bithigh_stdev,100)
toohigh_x = np.linspace(toohigh_mean-3*toohigh_stdev,toohigh_mean+3*toohigh_stdev,100)

plt.plot(toocheap_x, stats.norm.pdf(toocheap_x, toocheap_mean, toocheap_stdev),color='darkseagreen',lw=0.5)
plt.plot(goodbargain_x, stats.norm.pdf(goodbargain_x, goodbargain_mean, goodbargain_stdev),color='olivedrab',lw=0.5)
plt.plot(bithigh_x, stats.norm.pdf(bithigh_x, bithigh_mean, bithigh_stdev),color='lightsalmon',lw=0.5)
plt.plot(toohigh_x, stats.norm.pdf(toohigh_x, toohigh_mean, toohigh_stdev),color='orangered',lw=0.5)

# plt.axvline(x=toocheap_mean,color='k',lw=0.5)
# plt.axvline(x=goodbargain_mean,color='k',lw=0.5)
# plt.axvline(x=bithigh_mean,color='k',lw=0.5)
# plt.axvline(x=toohigh_mean,color='k',lw=0.5)
# plt.axvline(x=41,color='k',lw=0.5)
# plt.axvline(x=88,color='k',lw=0.5)
# plt.axvline(x=73,color='k',lw=0.5)

plt.plot([41,88],[0.0225,0.0225],color='navy',lw=0.85)
plt.plot([41,41],[0.0065,0.0225],color='k',lw=0.15)
plt.plot([73,73],[0.0116,0.0225],color='k',lw=0.15)
plt.plot([88,88],[0.0054,0.0225],color='k',lw=0.15)
plt.scatter([41,73,88],[0.023,0.023,0.023],marker='v')
plt.title('Product Pricing, 52 responses\nmin=41 avg=73 max=88')

plt.yticks([])
plt.xlabel("Price in USD")

plt.xlim(0,200)
plt.ylim(0,0.035)
plt.show()