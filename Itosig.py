
from olexFunctions import OlexFunctions
OV = OlexFunctions()

import os
import htmlTools
import olex
import olx
import gui
import time
import math
from PluginTools import PluginTools as PT
import sys
import plotly 
import plotly.graph_objects as go
from plotly.subplots import make_subplots as ms
import numpy as np


debug = bool(OV.GetParam("olex2.debug", False))
instance_path = OV.DataDir()

try:
  from_outside = False
  p_path = os.path.dirname(os.path.abspath(__file__))
except:
  from_outside = True
  p_path = os.path.dirname(os.path.abspath("__file__"))

z = open(os.sep.join([p_path, 'def.txt'])).readlines() 
d = {}
for line in z:
  line = line.strip()
  if not line or line.startswith("#"):
    continue
  d[line.split("=")[0].strip()] = line.split("=")[1].strip()

p_name = d['p_name']
p_htm = d['p_htm']
p_img = eval(d['p_img'])
p_scope = d['p_scope']

OV.SetVar('Itosig_plugin_path', p_path)

class Itosig(PT):

  def __init__(self):
    super(Itosig, self).__init__()
    self.p_name = p_name
    self.p_path = p_path
    self.p_scope = p_scope
    self.p_htm = p_htm
    self.p_img = p_img
    self.deal_with_phil(operation='read')
    self.print_version_date()
    if not from_outside:
      self.setup_gui()
    OV.registerFunction(self.itosigvsreflins,True,"Itosig")
    OV.registerFunction(self.test,True,"Itosig")
    # END Generated =======================================

  def itosigvsreflins(self):
    pass
  
  def test(self):
    """Plots various I to sigma and reflexion plots"""
    global hkl_complete_list
    hkl_complete_list = []
    a = open(OV.HKLSrc()).readlines()
    cell = [float(x) for x in olx.xf.au.GetCell().split(',')]
    
    if a[-1].startswith("E") == True :                 # Checking if INS header is attached, and removes ins header if so
      print("Ejecting INS header: last 16 lines")
      a = a[:len(a)-16]
    if a[-1].startswith("_") == True:
      print("Ejecting INS header: last 18 lines")
      a = a[:len(a)-18]    
    else:
      a = a[:len(a)-1] 

    hkl_comp1 = hkl_complete(a, cell)
    file_2 = olex.f('fileOpen("Choose comparison input HKL file", "HKL files|*.hkl",filepath())')
    if file_2 == "":                                  # Checking if a second file is loaded and proceeding with formatting
      pass
    else:
      b = open(file_2).readlines()
      if b[-1].startswith("E") == True:               # Checking if INS header is attached, and removes ins header if so
        b = b[:len(b)-16]
      if b[-1].startswith("_") == True:
        b = b[:len(b)-18]
      else:
        b = b[:len(b)-1]       
      print("END of hkl:", b[-1])
      hkl_comp2 = hkl_complete(b, cell)               
    
    sorted_inv_dspacings = []
    sorted_Itosig = []
    max_itosigs =  []
    nr_reflins = []
    for elem in hkl_complete_list:                    # generates a list of sorted d_spacings (lowest to highest)
      print(elem.hkl_comp)
      zet = ret_entrylist_of_listoflists(elem.hkl_comp, 6)
      sor_zet=[sorted(zet)]
      sorted_inv_dspacings.append(sor_zet)
      ito = ret_entrylist_of_listoflists(elem.hkl_comp, 5)
      sor_ito=[sorted(ito)]
      sorted_Itosig.append(sor_ito)
      max_itosigs.append(max_value(elem.hkl_comp, 5))
      nr_reflins.append(elem.max_values()[-1])
      print("Max Int, Max Sigma, Max I/sig, Avg_I/sig, total number of Reflexions used:", elem.max_values())

    try:
      plotting_everything(sorted_inv_dspacings, sorted_Itosig, 30, max(max_itosigs), nr_reflins)
    except:
      pass
            
def ret_entrylist_of_listoflists(listoflists, i):                       # help_funcs for further processing of the lists
    """Input: list of lists, output: list of float entries i of sublist"""
    ret = []
    for elem in listoflists:
      ret.append(float(elem[i]))
    return ret      
      
def max_value(inplist, i):
  """returns the max value of the i entry of a list"""
  return max([sublist[i] for sublist in inplist])  

def min_value(inplist, i):
  """returns the min value of the i entry of a list"""
  return min([sublist[i] for sublist in inplist])

def avg_list(inplist, i):
  """return the lists i entrys average"""
  counter = 0
  for sublist in inplist:
    counter += 1
  return (sum(sublist[i] for sublist in inplist))/counter

def calc_d_spacing(h, k, l, cell, Raumgruppe):
  """Calculated the respective d spacing. Params: h, k, l, cellparams, spacegroup"""
  a = cell[0]
  b = cell[1]
  c = cell[2]
  alphadeg = cell[3]
  betadeg = cell[4]
  gammadeg = cell[5]
  if h != 0 or k != 0 or l != 0:
    if Raumgruppe == "triklinic":
      alpha = alphadeg * (math.pi/180)
      beta = betadeg * (math.pi/180)
      gamma = gammadeg * (math.pi/180)
  
      h_a = h/a
      k_b = k/b
      l_c = l/c
  
      cos_alpha1 = math.cos(alpha)
      cos_beta1 = math.cos(beta)
      cos_gamma1 = math.cos(gamma)
  
      if abs(cos_alpha1) > 10**-6:
        cos_alpha = cos_alpha1
      else: 
        cos_alpha = 0
      if abs(cos_beta1) > 10**-6:
        cos_beta= cos_beta1
      else:
        cos_beta = 0
      if abs(cos_gamma1) > 10**-6:
        cos_gamma = cos_gamma1
      else:
        cos_gamma = 0
  
      A = [[h_a, cos_gamma, cos_beta], [k_b, 1, cos_alpha], [l_c, cos_alpha, 1]]
      B = [[1, h_a, cos_alpha], [cos_gamma, k_b, cos_alpha], [cos_beta, l_c, 1]]
      C = [[1, cos_gamma, h_a], [cos_gamma, 1, k_b], [cos_beta, cos_alpha, l_c]]
      D = [[1, cos_gamma, cos_beta], [cos_gamma, 1, cos_alpha], [cos_beta, cos_alpha, 1]]
  
      dstar = (h_a * la.det(A) + h_a * la.det(B) + l_c * la.det(C))*(1/la.det(C)) 
      dspace = (1/dstar)**0.5
      return dspace
    
    if Raumgruppe == "monoclinic":
      gamma = math.radians(betadeg)
      dh = (abs(h)**2/((a)**2*math.sin(gamma)**2))
      dk = (abs(k)**2)/((b)**2*math.sin(gamma)**2)-(2*(h)*(k)*math.cos(gamma))/(b**2*math.sin(gamma)**2)
      dl =  abs(l)**2/(c)**2
      dstar = dh + dk +dl
      dspace = (1/dstar)**0.5
      return dspace
    if Raumgruppe == "orthorhombic":
      dstar = h**2/a**2 + k**2/b**2 + l**2/c**2
      dspace = (1/dstar)**0.5
      return dspace
    if Raumgruppe == "tetragonal":
      dstar = (h**2 + k**2 + l**2)/((a/c)**2)/a**2
      dspace = (1/dstar)**0.5
      return dspace
    if Raumgruppe == "cubic":
      dstar = (h**2 + k**2 + l**2)/a**2
      dspace = (1/dstar)**0.5
      return dspace
    if Raumgruppe == "trigonal":
      print("this should be done by now, but isnt")
    if Raumgruppe == "hexagonal":
      dstar = (0.75*(h**2+k**2+h*k)+l**2*(a/c)**2)*1/a**2
      dspace = (1/dstar)**0.5
      return dspace

def gen_hkl_comp(inp, cell):
    """generates a list with the format of h, k, l, I, sig, I/sig, 1/d-spacing"""
    b = inp
    cell_params = cell
    SG = ""
    if cell[0] == cell[1] == cell[2] and cell[3] == cell[4] == cell[5] == 90.0:
      SG = "cubic"
    elif cell[0] == cell[1] == cell[2] and cell[3] == cell[4] == cell[5]:
      SG = "trigonal"
    elif cell[0] == cell[1] != cell[2] and cell[3] == cell[4] == 90.0 and cell[5] == 120.0:
      SG = "hexagonal"
    elif cell[0] == cell[1] != cell[2] and cell[3] == cell[4] == cell[5] == 90.0:
      SG = "tetragonal"
    elif cell[0] != cell[1] != cell[2] and cell[3] == cell[4] == cell[5] == 90.0:
      SG = "orthorhombic"
    elif cell[0] != cell[1] != cell[2] and cell[3] == cell[5] == 90.0 != cell[4]:
      SG = "monoclinic"
    else:
      SG = "triklinic"
    ret = []
    counter = 1
    for i in range(len(b)-1):
      if b[i][0] != "0" and b[i][1] != "0" and b[i][2] != "0" and b[i][3] != "0.00" and b[i][4] != "0.00" and b[i][5] != "0.00":
          c = []
          c = b[i][:4].split() + b[i][5:9].split() + b[i][9:13].split() + b[i][12:21].split() + b[i][21:29].split()
          c.append(float(b[i][12:21].split()[0]) / float(b[i][21:31].split()[0]))
          c.append(1/(calc_d_spacing(float(c[0]), float(c[1]), float(c[2]), cell, SG)))
          ret.append(c)
          counter += 1
      else: 
        continue       
    return [ret, counter]   

def plotting_everything(reflins_inv_dspac, itosig, bins, maxItosig, nrs):
  """Takes sorted lists of dspacing and I/sig values and creates plots of both vs numbers of reflexions"""
  fig = ms(rows = 2, cols = 2, subplot_titles=("d-spacing VS number of reflins", "Itosig vs Number of reflins", "Test3", "Test4"))
  count = 0
  
  for x in reflins_inv_dspac:
    if count == 0:
      name ="Alpha"
    else:
      name = "Beta"
    y = np.arange(len(x[0]))
    fig.add_trace(go.Scatter(x=x[0], y=y, name = name), row=1, col=1)
    count+=1
  
  count = 0 
  for elem in itosig:
    for x in reflins_inv_dspac:
      if count == 0:
        name ="Alpha"
      else:
        name = "Beta"
    rejections = 0
    for i in range(len(elem[0])):
      if elem[0][i] > 0.1:
        elem[0][i] = math.log10(elem[0][i])
      else:
        rejections += 1
    print('While plotting,', rejections, "reflexions have ben rejected due to very low I/sig (< 0.1) for the", name, "dataset")
        
    fig.add_trace(go.Histogram(x=elem[0], name = name, xbins= dict(start = 0, end = math.log10(maxItosig), size = math.log10(maxItosig)/bins)), row = 1, col = 2)
    count+=1
    
  fig.update_xaxes(title_text = "1/dspacing", row=1, col=1)
  fig.update_yaxes(title_text = "Number of reflexions", row=1, col=1)
  fig.update_xaxes(title_text = "binned log(I/sig)", row=1, col=2)
  fig.update_yaxes(title_text = "Number of reflexions", row=1, col=2)
  fig.update_layout(height=1500, width=2000)
  
  try:
    fig.show()
  except:
    print("plotting not successfull")

class hkl_complete:
       
    def __init__(self, file, cell_params):
      self.file = file
      self.cell_params = cell_params
      x = gen_hkl_comp(file, cell_params)
      self.hkl_comp = x[0]
      self.nr_reflins_used = x[1]
      self.hkl_comp.pop()
      self.max_int = max_value(self.hkl_comp, 3)
      self.max_sig = max_value(self.hkl_comp, 4)
      self.max_itosig = max_value(self.hkl_comp, 5)
      self.calc_avg_itosig = avg_list(self.hkl_comp, 6)
      hkl_complete_list.append(self)
     
    def max_values(self):
      return [self.max_int, self.max_sig, self.max_itosig, self.calc_avg_itosig, self.nr_reflins_used]
    
    def __call__(self):
      return self.hkl_comp

Itosig_instance = Itosig()

print("Ok")

