
**Please refer to "_ getting started.pdf" for instructions.**

This assignment comports different steps.

[1] Convert a Whitted-style ray tracer (No ray bounces / Reflection / Refraction & No acceleration structures), shooting primary rays and returning the color of the closest object into a Kajiya-style path tracer (With Reflection / Refraction and still no acceleration structures)
Step 1 commit: https://github.com/Fannollost/INFOMAGR_assignment_1/commit/555c9c2482aa66b78fa7b5315e557c354cff622c

Competence acquired :
* Shooting primary rays from the camera screen
* Rotate and translate the camera screen
* FOV change and fish-eye effect
* Basic primitives intersection (Planes, spheres, and triangles)
* Triangulated mesh load and intersection
* Intersection with a textured skydome
* Different light sources types (Spotlight, Point light, and area light)
* Different material types (Diffuse, dielectric, and reflective)
  
![Schermopname (34)](https://github.com/Fannollost/INFOMAGR_assignment_1/assets/47923208/4a078149-a946-488d-aeb2-70e8a75728ec)

[2] Accelerate a Kajiya-style path tracer (With Reflection / Refraction and still no acceleration structures) by building and using a BVH (Binding Volume Hierarchy) or a QBVH (QuadTree Binding Volume Hierarchy)

Competence acquired :
* Build, Refit, and intersect a BVH (Binding Volume Hierarchy)
* Different split heuristics (Longest axis median split, Longest axis middle split, Surface Area Heuristics, Binned Surface Area Heuristic)
* Analysis and comparison of different heuristics


![Average FPS](https://github.com/Fannollost/INFOMAGR_assignment_1/assets/47923208/7f9861e3-73ad-45f1-8cd0-87d5dcd0bd99)
![Treedepth](https://github.com/Fannollost/INFOMAGR_assignment_1/assets/47923208/545a6b5f-4821-4da3-aec9-e04b42f876f4)
![tri](https://github.com/Fannollost/INFOMAGR_assignment_1/assets/47923208/8aff6690-6e41-43cf-9daf-2224d3698545)


Step 2 final commit: https://github.com/Fannollost/INFOMAGR_assignment_1/commit/555c9c2482aa66b78fa7b5315e557c354cff622c

[3] Implement a recent paper: Here, use QLearning to influence the sampling direction (Dahm, K., & Keller, A. (2017). Learning Light Transport the Reinforced Way. CoRR, abs/1701.07403. http://arxiv.org/abs/1701.07403)

Competence acquired :
• Initialize sampling positions
• Pick sampling direction according to the QValue of neighboring points
• Store and update directions with a corresponding probability per sampling point

![debug](https://github.com/Fannollost/INFOMAGR_assignment_1/assets/47923208/d2121e6b-a842-4d0e-ad81-1630bd17735c)

Code is fully public domain. Use as you please.

The template author can be contacted at bikker.j@gmail.com.

