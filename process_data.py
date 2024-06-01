import csv
import matplotlib.pyplot as plt
import numpy as np

def pwr_sum(readings):
	sum = 0
	for i in range(10, len(readings)-10):
		sum += readings[i]
	return sum

cores = 1

while cores <= 64:
	x = 4
	y = 4
	z = 4
	while x*y*z <= 64*64*64:
		case = str(cores) + '-' + str(x) + '-' + str(y)+ '-' + str(z)

		total = []
		core0 = []
		mem0 = []
		sram0 = []
		soc0 = []
		core1 = []
		mem1 = []
		sram1 = []
		soc1 = []
		
		with open('hpcg-' + case + '/pwr-' + case + '.csv', newline='') as csvfile:
			reader = csv.DictReader(csvfile)
			for row in reader:
					core0.append(float(row['pwr_core0']))
					mem0.append(float(row['pwr_mem0']))
					sram0.append(float(row['pwr_sram0']))
					soc0.append(float(row['pwr_soc0']))
					core1.append(float(row['pwr_core1']))
					mem1.append(float(row['pwr_mem1']))
					sram1.append(float(row['pwr_sram1']))
					soc1.append(float(row['pwr_soc1']))
					total.append(float(row['pwr_core0']) + float(row['pwr_mem0']) + float(row['pwr_sram0']) + float(row['pwr_soc0']) + float(row['pwr_core1']) + float(row['pwr_mem1']) + float(row['pwr_sram1']) + float(row['pwr_soc1']))
		
		lines = []
		speed = 0.0
		memUsed = 0.0
		flops = 0

		with open('hpcg-' + case + '/HPCG-Benchmark.txt') as f:
			lines = f.readlines()

		for line in lines:
			if 'result is VALID' in line:
				speed = float(line.rpartition('=')[2])
				# print(speed)
			if 'Floating Point Operations Summary::Total' in line:
				flops = float(line.rpartition('=')[2])
				# print(speed)
			if 'Total memory used' in line:
				memUsed = float(line.rpartition('=')[2])
				# print(memUsed)

		with open('results.csv', 'a', newline='') as csvfile:
			fnames = ['Test Case', 'core0', 'mem0', 'sram0', 'soc0', 'core1', 'mem1', 'sram1', 'soc1', 'total', 'speed', 'memUsed', 'flops']
			writer = csv.DictWriter(csvfile, fieldnames=fnames)
			# print(speed)
			if(case == '1-4-4-4'):
				writer.writeheader()
			writer.writerow({'Test Case' : case, 'core0' : str(pwr_sum(core0)), 'mem0' : str(pwr_sum(mem0)), 'sram0' : str(pwr_sum(sram0)), 'soc0' : str(pwr_sum(soc0)), 'core1' : str(pwr_sum(core1)), 'mem1' : str(pwr_sum(mem1)), 'sram1' : str(pwr_sum(sram1)), 'soc1' : str(pwr_sum(soc1)), 'total' : str(pwr_sum(total)), 'speed' : speed, 'memUsed': memUsed, 'flops' : flops})	

		if z < y:
			z = z*2
		elif y < x:
			y = y*2
		else:
			x = x*2
	cores = cores*2

total = np.array(total)


xpoints = np.array(np.arange(start=0, stop=len(total), step=1))
ypoints = total
plt.plot(xpoints, total)
# plt.show()
