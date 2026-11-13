#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:34:19 2021

@author: kendrick
"""

import numpy as np

# compute unknown displacements 
def ComputeDisplacements(K, F, n_unknowns):
    # extract submatrix of unknowns
    K11 = K[0:n_unknowns,0:n_unknowns]
    F1 = F[0:n_unknowns]
    
    d = np.linalg.solve(K11,F1)
    
    return d

# postprocess the forces at known displacement nodes
def PostprocessReactions(K, d, F, n_unknowns, nodes):
    # These are computed net forces and do not
    # take into account external loads applied
    # at these nodes
    F = np.matmul(K[n_unknowns:,0:n_unknowns], d)
    
    # Postprocess the reactions
    for node in nodes:
        if node.xidx >= n_unknowns:
            node.AddReactionXForce(F[node.xidx-n_unknowns][0] - node.xforce_external)
        if node.yidx >= n_unknowns:
            node.AddReactionYForce(F[node.yidx-n_unknowns][0] - node.yforce_external)
        
    return F

# determine internal member loads
def ComputeMemberForces(bars):
    for bar in bars:
        E = bar.E
        A = bar.A
        L = bar.Length()
        lam_x, lam_y = bar.LambdaTerms()

        
        near_node = bar.init_node
        far_node = bar.end_node
        
        d_nx = near_node.xdisp
        d_ny = near_node.ydisp
        d_fx = far_node.xdisp
        d_fy = far_node.ydisp
        
        array_1 = np.array([-lam_x, -lam_y, lam_x, lam_y])
        array_2 = np.array([d_nx, d_ny, d_fx, d_fy])
        q_F = (A * E/L) * array_1 @ array_2
        bar.axial_load = q_F
    # COMPLETE THIS FUNCTION
    # Compute member forces for all bars using equation 14-23 
    return bars
    
# compute the normal stresses
def ComputeNormalStresses(bars):
    # COMPLETE THIS FUNCTION
    for bar in bars:
        if not np.isnan(bar.axial_load):
            bar.normal_stress = bar.axial_load / bar.A
    # Compute normal stress for all bars
    return bars

# compute the critical buckling load of a member
def ComputeBucklingLoad(bars):
    # COMPLETE THIS FUNCTION
    for bar in bars:
        E = bar.E
        I = bar.GetWeakSecondMoment()
        L =bar.Length() *12
        K = 1
        
        if E > 0 and I > 0 and L > 0:
            P_cr = (np.pi**2 * E * I) / ((K * L)**2)
            bar.buckling_load = P_cr
    # Compute critical buckling load for all bars
    return bars
