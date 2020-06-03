
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
import matplotlib as plt

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
    
    if a[-1].startswith("E") == True:                 # Checking if INS header is attached, and removes ins header if so
      a = a[:len(a)-12]
      
    cell = [float(x) for x in olx.xf.au.GetCell().split(',')]
    hkl_comp1 = hkl_complete(a, cell)
    file_2 = olex.f('fileOpen("Choose input HKL file", "HKL files|*.hkl",filepath())')
    
    if file_2 == "":                                  # Checking if a second file is loaded and proceeding with formatting
      pass
    else:
      b = open(file_2).readlines()
      if b[-1].startswith("E") == True:               # Checking if INS header is attached, and removes ins header if so
        b = b[:len(b)-12]      
      hkl_comp2 = hkl_complete(b, cell)               # For now, the compare file uses the EXACT same cell paramenters as the first file
      
    for elem in hkl_complete_list:                    # generates a list of sorted d_spacings (lowest to highest)
      zet = ret_entrylist_of_listoflists(elem.hkl_comp, 6)
      sor_zet=[sorted(zet)]
      print(sor_zet)
      
def ret_entrylist_of_listoflists(listoflists, i):
    """Input: list of lists, output: list of entries i of sublist"""
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
  return (sum(float(sublist[i]) for sublist in inplist))/counter

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
      print("this should be done by now, but isnt")    

def gen_hkl_comp(inp, cell):
    """generates a list with the format of h, k, l, I, sig, I/sig, d-spacing"""
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
    for i in range(len(b)-1):
      if b[i][0] != "0" and b[i][1] != "0" and b[i][2] != "0" and b[i][3] != "0.00" and b[i][4] != "0.00" and b[i][5] != "0.00":
        c = []
        c = b[i][:4].split() + b[i][5:9].split() + b[i][9:12].split() + b[i][12:21].split() + b[i][21:30].split()
        
        c.append(float(b[i][12:21].split()[0]) / float(b[i][21:31].split()[0]))
        c.append(calc_d_spacing(float(c[0]), float(c[1]), float(c[2]), cell, SG))
        ret.append(c)             
      else: 
        continue       
    return ret   

  

      

class hkl_complete:
       
    def __init__(self, file, cell_params):
      self.file = file
      self.cell_params = cell_params
      self.hkl_comp = gen_hkl_comp(file, cell_params)
      self.hkl_comp.pop()
      self.max_int = max_value(self.hkl_comp, 3)
      self.max_sig = max_value(self.hkl_comp, 4)
      self.max_itosig = max_value(self.hkl_comp, 5)
      self.calc_itosig = avg_list(self.hkl_comp, 6)
      hkl_complete_list.append(self)
    
    
    def reflIns_vs_dspacing_plot(self):
      uns = []
      uns.append(ret_entrylist_of_listoflists(self.hkl_comp[6]))
      print(uns)
      return
      

Itosig_instance = Itosig()

print("OK.")

