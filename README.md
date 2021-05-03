# AeroDesign-Performance

Performance study of a small aircraft, builded to compete in the SAE Brasil AeroDesign 2018.

The whole aircraft design process was aided by python and fortran scripts, and validate on real flights.

The **performance.py** script aim to test the performance of a concptual aircraft design. For this task the algorithm propused by *Roskam* was used to calculate the theorical drag coefficient(Cd), and plot a *V-N Diagram*.

### Long Story Short 

The script generate the curve a O.S 55 engine measured in wind tunnel and compare with the *Roskam algorithm* results, as the **Figure 01** shows, it's possbile to find out the maximum and de minimum flight speed. The **Figure 02** ilustrate the resulting *V-N Diagram*. The final results are showed next the **Figure 02**.

![performance](/images/thurst_performance.png)
**Figure 01**
![v_n](/images/v_n.png)
**Figure 02**
```
---------------- Takeoff Distance --------------------

Takeoff Distance:  123.76 m

Maximum Climb Angle:  4.72 °

CL/CD Ratio  12.12

Maximum Climb Ratio:  3.0 m/s

----------------------- Score -----------------

Empty Mass:  2.5 kg

Load:  6.0 kg

Total Mass  8.5 kg

Strucutral Efficiency:  2.4

Score:  42.0 Points
---------------- Cruise Flight --------------------

Maximum Velocity:  39.4 m/s error +-  5 %

Diving Velocity:  32.5 m/s
```
