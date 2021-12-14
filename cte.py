# -*- coding: utf-8 -*-
"""
Created on Sun May 23 16:50:27 2021

@author: rodri
"""

import numpy as np

Rv=461.51
cp=1004
T_lt = 273
e_lt = 611.2
Rd=287.053
epsi = Rd/Rv
p_min= 1.03e5
p_top = 1e4
p = np.linspace(p_min,p_top,int((p_min - p_top)/500 + 1 ))
T_min = -150
T_max = 40
N = (T_max - T_min)//1 +1

T=np.linspace(T_min,T_max,N)
beta = -40

ws_color = '#0d3c69'