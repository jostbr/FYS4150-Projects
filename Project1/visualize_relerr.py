
import numpy as np
import matplotlib.pyplot as plt

problem_sizes = np.logspace(1, 7, 7, endpoint = True).astype(int)
step_sizes = 1.0/(problem_sizes + 1)

max_err_list = list()

for i in range(len(problem_sizes)):
    max_err = 0.0

    with open("../build-Project1/result{}.txt".format(problem_sizes[i]), "r") as result_file:
        for line_num, line_string in enumerate(result_file):
            if (line_num > 1 and line_num < problem_sizes[i] + 2):
                words = line_string.split()
                rel_err = abs(float(words[-1]))

                if (rel_err > max_err):
                    max_err = rel_err

        max_err_list.append(max_err)

print(problem_sizes)
print(np.log10(np.array(max_err_list)))

plt.figure(figsize = (8, 7))
plt.style.use("seaborn")
plt.plot(problem_sizes, max_err_list, linewidth = 2)
plt.loglog()
plt.xlabel("Logarithmic number of grid points, N", fontname = "serif", fontsize = 13)
plt.ylabel("Logarithmic Relative error [dim-less]", fontname = "serif", fontsize = 13)
plt.title("Relative error from $N=10$ to $N=10^7$", fontname = "serif", fontsize = 17)
plt.show()



