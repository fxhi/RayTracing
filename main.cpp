// Header OpenMP
#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

// Header C++
#include <iostream>
#include <fstream>

// Header Ray tracing
#include "rtweekend.h"

#include "camera.h"
#include "bvh.h"
#include "color.h"
// #include "hittable_list.h"
#include "material.h"
#include "sphere.h"
#include "moving_sphere.h"



color ray_color(const ray& r, const hittable& world, int depth)
{
    hit_record rec;

    // If we've exceed the ray bounce limit, no more light is gathered.
    if ( depth <= 0)
        return color(0,0,0);

    if(world.hit(r, 0.001, infinity, rec)) // 0.001 remove shadow acne ( and faster the code !)
    {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r,rec,attenuation,scattered))
            return attenuation * ray_color(scattered, world, depth-1);
        return color(0,0,0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}

hittable_list random_scene() {
    hittable_list world;

    hittable_list objects;

    auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0,-1000,0), 1000, ground_material));

    for (int a = -11; a < 11; a++)
    {
        for (int b = -11; b < 11; b++)
        {
            auto choose_mat = random_double();
            point3 center( a + 0.9*(fabs(a)/10.0 - 0.1), 0.2, b + (fabs(b)/10.0 - 0.1));
            // point3 center( a + 0.9*random_double(), 0.2, b + 0.9*random_double());
            
            // fix value between 0 and 1
            auto val1 = fabs((fabs(a)/10.0 - 0.1)-0.5);
            auto val2 = (fabs(b)/10.0 - 0.1);
            auto val3 = fabs((fabs(b)/10.0 - 0.1)-0.5);

            if ((center - point3(4, 0.2, 0)).length() > 0.9)
            {
                shared_ptr<material> sphere_material;

                if(choose_mat < 0.8)
                {
                    //diffuse
                    auto albedo = color(val1, val2, val3) * color(val1, val2, val3);
                    // auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    auto center2 = center + vec3(0, random_double(0, 0.5), 0);
                    world.add(make_shared<moving_sphere>(
                        // center, center2, 0.0, 0.1, 0.2, sphere_material));
                        center, center2, 0.0, 1.0, 0.2, sphere_material));
                    // world.add(make_shared<sphere>(
                    //     center,  0.2, sphere_material));
                }
                else if ( choose_mat < 0.95) 
                {
                    auto albedo = color(val1, val2, val3)*0.5 + color(0.5, 0.5, 0.5);
                    auto fuzz = val1*0.5;
                    // auto albedo = color::random(0.5, 1);
                    // auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else
                {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.5);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    objects.add(make_shared<bvh_node>(world, 0, 1));

    return objects;
    // return world;
}



int main() {

    //> --- Redirecting cout flux to a file

    std::ofstream out_file{"imageRayTracingOut.ppm"};
    std::streambuf *coutbuf {std::cout.rdbuf()}; // save old buff
    std::cout.rdbuf(out_file.rdbuf()); // redirectng std::cout to out variable

    //>--- Image parameters

    // const auto aspect_ratio = 4.0 / 3.0;
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 50;
    const int max_depth = 8;

    //> --- World

    auto world = random_scene();

    //> --- Camera

    point3 lookfrom(13,2,3);
    point3 lookat(0,0,0);
    vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;

    // camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);
    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);


    //> --- Render

    // Image are stored with the PPM format

    // Header of the ppm file
    std::cout << "P3\n" << image_width << ' ' << image_height << ' ' << "\n255\n";

    // Draw image pixels from top to bottom, left to right

    //For OpenMP parallelization
    int array_ir[image_width][image_height], array_ig[image_width][image_height] , array_ib[image_width][image_height];
    int ir, ig, ib;

    // int n;
    // #pragma omp parallel for
    // for(n=0;n<8;n++){
    //     printf("Element %d traité par le thread %d \n",n,omp_get_thread_num());
    // }
    
    // #pragma omp parallel for
    std::cerr << "Number of threads set for OpenMP wich can be used" << omp_get_max_threads() << std::endl;
    int j;
    int k = 0;
    #pragma omp parallel for num_threads(4) private(ir, ig, ib) shared(j,array_ir, array_ig, array_ib)
    for (j = image_height-1; j >= 0; j--) {
    // for (int j = image_height-1; j >= 0; j--) {
        #pragma omp critical
        {
            std::cerr << "\rNumber of row treated : " << k << "/" << image_height  << " thread " << omp_get_thread_num() << std::flush;
        }
        // std::cerr << "\rRemaning rows : " << j+1 << "/" << image_height  << " thread " << omp_get_thread_num() << std::endl;
        // std::cerr << "\rRemaning rows : " << j+1 << "/" << image_height  << ' ' << std::flush;
        // std::cerr << "\rScanlines remaining: " << j << ' ' << std::endl;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; s++)
            {
                auto u = (i + random_double()) / (image_width-1);
                auto v = (j + random_double()) / (image_height-1);

                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(ir, ig, ib, pixel_color, samples_per_pixel);
            array_ir[i][j] = ir;
            array_ig[i][j] = ig;
            array_ib[i][j] = ib;
            // write_color(std::cout, pixel_color, samples_per_pixel);
        }
        k++;
    }

    std::cerr << std::endl;
    std::cerr << "Writing file..." << std::endl;
    for (int j = image_height-1; j >= 0; j--) {
        for (int i = 0; i < image_width; i++)
        {
            std::cout << array_ir[i][j] << ' ' << array_ig[i][j] << ' ' << array_ib[i][j] <<  '\n';
        }
    }
    std::cerr << "Done." <<  std::endl;

    //> --- Restoring cout flux

    std::cout.rdbuf(coutbuf);        // restore cout original streambuf
    out_file.close();

    return 0;
}
