import numpy as np
import math
import seqlogo
import matplotlib.pyplot as plt
import os
import sys


current_script_path = sys.argv[0]
subpackage_folder_path = os.path.dirname( __file__ )
if subpackage_folder_path not in sys.path:
    sys.path.append(subpackage_folder_path)

sys.path.append("/Users/student/Downloads/weblogo-master")
import weblogo as wl

import glob_var

import modify_seed
import type_conversions


def generate_frequencies(inp_w_motif):
    frequencies = np.zeros((inp_w_motif.linear_length, 4))
    inp_w_motif.get_linear_sequence()
    for i, nt in enumerate(inp_w_motif.linear_sequence):
        set_of_letters = glob_var._degenerate_nts_mapping[nt]
        for let in set_of_letters:
            column_index = glob_var._rna_alphabet_nt_mapping[let]
            frequencies[i, column_index] = 1
        if len(set_of_letters) == 0:
            frequencies[i, :] = np.array([1,1,1,1])

    return frequencies


def process_freq_for_seqlogo(frequencies):
    background_NA_dict = {nt: 0.25 for nt in glob_var._rna_alphabet_list}
    background_NA_list = np.array(list(background_NA_dict.values()))
    # create PFM object for seqlogo
    pfm = seqlogo.Pfm(frequencies, alphabet_type='RNA', background=background_NA_list)
    # compute missing data: PPM, PWM, information content etc
    # no need to do weird reformatting between PFM, PPM, PWM
    # see https://pypi.org/project/seqlogo/
    # don't look at "Generate some frequency data and convert to PWM" - it's misleading
    PM = seqlogo.CompletePm(pfm, alphabet=glob_var._rna_alphabet_string)

    return PM


def draw_weblogo(inp_w_motif, out_file):
    frequencies = generate_frequencies(inp_w_motif)
    plottable_matrix = process_freq_for_seqlogo(frequencies)

    seqlogo.seqlogo(plottable_matrix, color_scheme='classic',
                    format='eps', size='medium', filename=out_file)


# def print_motif(inp_w_motif, do_print = True):
#     if do_print:
#         print("Short motif representation:")
#         inp_w_motif.print()
#         print("Full motif representation:")
#         inp_w_motif.print_linear()
#         print()



def playground(w_motifs_list):
    fig, axes = plt.subplots(ncols=2, nrows=len(w_motifs_list))

    for i, axes_row in enumerate(axes):
        for k, ind_axes in enumerate(axes_row):
            ind_axes.set_title('Axes %d %d' % (i, k))
            ind_axes.axis('off')
            draw_one_motif_v2(ind_axes, w_motifs_list[i])
            #ind_axes.plot(w_motifs_list[i].sequence)

    # for i in range(len(w_motifs_list)):
    #     current_ax = axes[i]
    #     current_ax.plot(w_motifs_list[i].sequence)

    plt.show()


def draw_secondary_structure_dev(w_motifs_list):
    cnt = len(w_motifs_list)
    maxlen = 0

    for motif in w_motifs_list:
        motif.adjust_linear_length()
        if motif.linear_length > maxlen:
            maxlen = motif.linear_length

    col = 5
    maxlen += 5
    xbase = 140
    ybase = 120
    xsize = xbase * 3 +maxlen / 3 * 30 * col
    ysize = ybase * 2 +maxlen / 2 * 30 * cnt /col

    print(xsize, ysize)

    #drawRNA(current_ax, xbase, ybase, structure, bases, l, r, lx, ly, rx, ry, colors={})



# def draw_one_motif(current_ax, w_motif):
#     sequence_string = w_motif.print_linear_sequence(return_string=True)
#     structure_string = w_motif.print_linear_structure(return_string = True)
#     print(structure_string)
#
#     xbase = 1
#     ybase = 1
#
#     drawRNA(current_ax, structure_string, sequence_string, l=0, r=len(structure_string), lx=0, ly=0, rx=1, ry=0, colors={})
#
#     # l = 0
#     # r = w_motif.linear_length
#     #
#     # L = 0
#     # count = 2
#     # for i in range(l+1, r):
#     #     if structure_string[i] == ">":
#     #         L -= 1
#     #     if L == 0:
#     #         count += 1
#     #     if structure_string[i] == "<":
#     #         L += 1
#     #
#     # print(count, L)
#
#
#
#
# def drawRNA(current_ax, structure, bases, l, r, lx, ly, rx, ry, colors):
#     print(l, r, lx, ly, rx, ry)
#
#     L = 0
#     count = 2
#     for i in range(l+1, r):
#         if structure[i] == ">":
#             L -= 1
#         if L == 0:
#             count += 1
#         if structure[i] == "<":
#             L += 1
#
#     th = 2 * 3.14159 / count
#     R = 1 / (2*math.sin(th/2))
#     h = R * math.cos(th/2)
#     (cx, cy) = (((lx+rx)/2.0)+h*(ly-ry), ((ly+ry)/2.0)+h*(rx-lx))
#     deg = math.atan2(ly-cy,lx-cx)
#
#
#     (i2,x2,y2) = (l,lx,ly)
#
#     for i in range(l+1, r):
#         if structure[i] == ">": # TODO: something about this is wrong
#             L -= 1
#         if L == 0:
#             deg -= th
#             (x,y) = (cx+R*math.cos(deg), cy+R*math.sin(deg))
#             # p.setcolour(118, 118, 118)
#             # p.setlinewidth (0.5)
#
#             print("X: ", deg, round(x,2), round(x2,2))
#             current_ax.plot([x2 * 10, x * 10], [y2 * 10, y * 10])
#
#
#             # p.setlinewidth (0.1)
#             # p.setcolour("black")
#             if structure[i] == ">":
#                 drawRNA(current_ax, structure, bases, i2, i, x2, y2, x, y, colors)
#             # col = colors->{bases->[i2]}
#             # p.setcolour(@col)
#             # p->circle({filled=>1},x2*10+xbase,y2*10+ybase,3)
#             # p->setcolour("black")
#             # p->circle({filled=>0},x2*10+xbase,y2*10+ybase,3)
#             # p->text({align=>"center", rotate=>"0"}, x2*10+xbase,y2*10+ybase-1.5, bases->[i2])
#             (x2,y2)=(x,y)
#             i2 = i
#         if structure[i] == "<":
#             L += 1
#     # p.setcolour(118, 118, 118)
#     # p.setlinewidth (0.5)
#     # current_ax.plot([x2 * 10, x * 10], [y2 * 10, y * 10])
#     # p.line(x2*10+xbase,y2*10+ybase,rx*10+xbase,ry*10+ybase)
#     # p.setlinewidth (0.1)
#     # p.setcolour("black")
#     # #my col = colors->{bases->[r-1]}
#     # #p->setcolour(@col)
#     # p.circle({filled=>1},x2*10+xbase,y2*10+ybase,3)
#     # p.setcolour("black")
#     # p.circle({filled=>0},x2*10+xbase,y2*10+ybase,3)
#     # p.text({align=>"center", rotate=>"0"}, x2*10+xbase,y2*10+ybase-1.5, bases->[r-1])
#     # p.setcolour(118, 118, 118)
#     # p.line(lx*10+xbase,ly*10+ybase,rx*10+xbase,ry*10+ybase) unless (lx==0 || ly==0)   ## 3-5 pair
#     # p.setcolour("black")
#     # my col = colors->{bases->[r]}
#     # p->setcolour(@col)
#     # #p->box({filled=>1},rx*10-2.5+xbase,ry*10+ybase-2,rx*10+2.5+xbase,ry*10+ybase+3)
#     # p->circle({filled=>1},rx*10+xbase,ry*10+ybase,3)
#     # p->setcolour("black")
#     # #p->box({filled=>0},rx*10-2.5+xbase,ry*10+ybase-2,rx*10+2.5+xbase,ry*10+ybase+3)
#     # p->circle({filled=>0},rx*10+xbase,ry*10+ybase,3)
#     #
#     # p->text({align=>"center", rotate=>"0"}, rx*10+xbase,ry*10+ybase-1.5, bases->[r])
#     # #p->text({align=>"center", rotate=>"90"}, rx*10+xbase+1.5,ry*10+ybase+0.5, bases->[r])
#     # }


def draw_one_motif_v2(current_ax, w_motif):
    sequence_string = w_motif.print_linear_sequence(return_string=True)
    structure_string = w_motif.print_linear_structure(return_string = True)
    print(structure_string)

    xbase = 1
    ybase = 1

    drawRNA_v2(current_ax, structure_string, sequence_string, l=0, r=len(structure_string), lx=0, ly=0, rx=1, ry=0, colors={})

    # l = 0
    # r = w_motif.linear_length
    #
    # L = 0
    # count = 2
    # for i in range(l+1, r):
    #     if structure_string[i] == ">":
    #         L -= 1
    #     if L == 0:
    #         count += 1
    #     if structure_string[i] == "<":
    #         L += 1
    #
    # print(count, L)




def drawRNA_v2(current_ax, structure, bases, l, r, lx, ly, rx, ry, colors):
    print(l, r, lx, ly, rx, ry)

    L = 0
    count = 2
    for i in range(l+1, r):
        print("Intermediate: %d, %d" % (count, L))
        if structure[i] == ">":
            L -= 1
        if L == 0:
            count += 1
        if structure[i] == "<":
            L += 1

    print("Count: %d, L: %d" % (count, L))
    print()

    th = 2 * 3.14159 / count
    R = 1 / (2*math.sin(th/2))
    h = R * math.cos(th/2)
    (cx, cy) = (((lx+rx)/2.0)+h*(ly-ry), ((ly+ry)/2.0)+h*(rx-lx))
    deg = math.atan2(ly-cy,lx-cx)


    (i2,x2,y2) = (l,lx,ly)

    for i in range(l+1, r):
        if structure[i] == ">": # TODO: something about this is wrong
            L -= 1
        if L == 0:
            deg -= th
            (x,y) = (cx+R*math.cos(deg), cy+R*math.sin(deg))
            # p.setcolour(118, 118, 118)
            # p.setlinewidth (0.5)

            current_ax.plot([x2 * 10, x * 10], [y2 * 10, y * 10])


            # p.setlinewidth (0.1)
            # p.setcolour("black")
            if structure[i] == ">":
                drawRNA_v2(current_ax, structure, bases, i2, i, x2, y2, x, y, colors)
            # col = colors->{bases->[i2]}
            # p.setcolour(@col)
            # p->circle({filled=>1},x2*10+xbase,y2*10+ybase,3)
            # p->setcolour("black")
            # p->circle({filled=>0},x2*10+xbase,y2*10+ybase,3)
            # p->text({align=>"center", rotate=>"0"}, x2*10+xbase,y2*10+ybase-1.5, bases->[i2])
            (x2,y2)=(x,y)
            i2 = i
        if structure[i] == "<":
            L += 1
