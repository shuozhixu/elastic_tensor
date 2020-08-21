The program, `elastic_tensor.c`, transforms a stiffness tensor for a given coordinate system to a stiffness tensor and a compliance tensor for another coordinate system. It requires an input file, `elas.in`, in the same directory.

The first line of `elas.in` is an integer, which can be 0, 1, 2, or 3. When it is 0, elastic anisotropy is assumed. When it is 1, elastic isotropy in Voigt form is assumed. When it is 2, elastic isotropy in Reuss form. When it is 3, elastic isotropy in Hill form is assumed.

Each of the next three lines of elas.in contains three real numbers. The first, second, and third line, respectively, corresponds to the crystallographic orientations of the _x_, _y_, and _z_ axis of the coordinate system based on which the input stiffness tensor is calculated, e.g.,

	1. 1. -2.
	-1. 1. 0.
	1. 1. 1.

The next six lines of `elas.in` contain 21 real numbers which are components of the input stiffness tensor, i.e,

<pre>
	<i>C</i><sub>11</sub>
	<i>C</i><sub>12</sub> <i>C</i><sub>22</sub>
	<i>C</i><sub>13</sub> <i>C</i><sub>23</sub> <i>C</i><sub>33</sub>
	<i>C</i><sub>14</sub> <i>C</i><sub>24</sub> <i>C</i><sub>34</sub> <i>C</i><sub>44</sub>
	<i>C</i><sub>15</sub> <i>C</i><sub>25</sub> <i>C</i><sub>35</sub> <i>C</i><sub>45</sub> <i>C</i><sub>55</sub>
	<i>C</i><sub>16</sub> <i>C</i><sub>26</sub> <i>C</i><sub>36</sub> <i>C</i><sub>46</sub> <i>C</i><sub>56</sub> <i>C</i><sub>66</sub>
</pre>

Each of the next three lines of `elas.in` contains three real numbers. The first, second, and third line, respectively, corresponds to the crystallographic orientations of the _x_, _y_, and _z_ axis of the coordinate system based on which the output stiffness and compliance tensors are calculated, e.g.,

	1. 0. 0.
	0. 1. 0.
	0. 0. 1.

To compile the program:

	gcc elastic_tensor.c -o elastic_tensor -lm

Then to run the program:

	./elastic_tensor

which creates an output file, `elas.out`. The file `elas.out` contains a stiffness tensor, which is Equation 9 of [a paper in _Intermetallics_](http://dx.doi.org/10.1016/j.intermet.2020.106844), and a compliance tensor.

If you use this program in your published work, please cite:

Shuozhi Xu, Yanqing Su, Irene J. Beyerlein, [Modeling dislocations with arbitrary character angle in face-centered cubic transition metals using the phase-field dislocation dynamics method with full anisotropic elasticity](http://dx.doi.org/10.1016/j.mechmat.2019.103200), Mech. Mater. 139 (2019) 103200
