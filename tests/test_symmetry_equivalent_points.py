# should match `plot_symmetry_equivalent_points` test in `symmetry.rs`
from pymatgen.symmetry.groups import SpaceGroup
import json
from matplotlib import pyplot as plt

point = [0.0, 0.0, 0.0]

x = []
y = []
z = []

for s in SpaceGroup.from_int_number(227).symmetry_ops:
    new_point = s.operate(point)
    x.append(new_point[0])
    y.append(new_point[1])
    z.append(new_point[2])

# remove duplicates
xyz = list(zip(x, y, z))
xyz = list(set(xyz))
x, y, z = zip(*xyz)

print(x.__len__())

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(x, y, z)
plt.show()
