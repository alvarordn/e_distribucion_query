# Libraries
import numpy as np

# Classes for grid generation
class grid:
    def __init__(self, nodes, lines, meas):
        self.nodes = self.add_nodes(nodes)                                      
        self.lines = self.add_lines(lines, self.nodes)   
        self.n = len(self.nodes)*2 - 1
        self.meas = self.add_meas(meas, self.nodes, self.lines, self.n)
        self.H = np.zeros((len(self.meas), self.n))
                
    def add_nodes(self, nodes):
        nodes_list = list()
        for item in nodes:
            nodes_list.append(node(item['id']))
        nodes_list[0].pointer = [None, len(nodes_list) - 1]
        for index, item in enumerate(nodes_list[1:]):
            item.pointer = [index, len(nodes_list) + index]
        return nodes_list
        
    def add_lines(self, lines, nodes):
        lines_list = list()
        for item in lines:
            lines_list.append(line(item['id'], item['From'], item['To'], item['R'], item['X'], nodes))        
        return lines_list
    
    def add_meas(self, meas, nodes, lines, n):
        meas_list = list()
        for item in meas:
            meas_list.append(measurement(item['id'], item['node'], item['line'], item['type'], item['value'], item['std'], nodes, lines, n))
        return meas_list
    
    def state_estimation(self, x, tol = 1e-5, niter = 100):
        x = np.array(x)
        self.build_W()
        self.assign(x)
        self.build_H(x)
        self.build_G()
        res, sol = list(), list()                
        self.compute_res(x)        
        res.append(self.res)   
        sol.append(x)
        iteration = 1
        print('')
        while (np.max(np.abs(res[-1])) > tol) and (iteration < niter):
            print(f'Iteration {iteration}, residual {np.max(np.abs(res[-1])):.3f}')
            x = self.update_x(x)
            sol.append(x)             
            self.compute_res(x)        
            res.append(self.res)
            self.build_H(x)
            self.build_G()
            iteration += 1
        self.compute_mags()
        print('')
        return res, sol
        
    def compute_mags(self):
        for node in self.nodes:
            node.Vx = complex(node.V*np.cos(node.theta), node.V*np.sin(node.theta))
        for line in self.lines:
            V1 = line.nodes[0].Vx
            V2 = line.nodes[1].Vx
            line.Ix = ((V1 - V2) / line.Z)
            line.I = np.abs(line.Ix)
        for node in self.nodes:
            node.Ix = np.sum([line.Ix if node == line.nodes[0] else -line.Ix for line in node.lines])
            node.I = np.abs(node.Ix)        
        
    def build_H(self, x):
        self.assign(x)
        for index, m in enumerate(self.meas):
            m.compute_h()
        self.H = np.array([item.dh for item in self.meas])
        
    def build_W(self):
        self.W = np.diag([1/(item.std**2) for item in self.meas])
        self.R = np.linalg.inv(self.W)
        
    def build_G(self):
        self.G = self.H.T.dot(self.W).dot(self.H)
     
    def assign(self, x):
        for index in range(len(self.nodes) - 1):
            self.nodes[1 + index].theta = x[index]        
        for index in range(len(self.nodes)):
            self.nodes[index].V = x[index + len(self.nodes) - 1]
            
    def compute_res(self, x):
        self.assign(x)
        self.res = [m.value - m.h() for m in self.meas]        
        
    def norm_res(self):
        M = self.R - self.H.dot(np.linalg.inv(self.G).dot(self.H.T))
        self.res_norm = [item[0]/np.sqrt(item[1]) for item in zip(np.abs(self.res), np.diag(np.abs(M)))]
        
    def update_x(self, x):
        return x + np.linalg.solve(self.G, self.H.T.dot(self.W).dot(self.res)) 
    
    def report(self):
        for node in self.nodes:
            print(f'Node {node.ref}: U = {node.V:.3f},\t  theta = {node.theta*180/np.pi:.3f}, \t Ix = {node.Ix:.3f}')
        print('')
        for line in self.lines:
            print(f'Line {str(line.nodes[0].ref)+"-"+str(line.nodes[1].ref)}: Ix = {line.Ix:.3f}')
        print('')
        print('\t\t meas. \t\t state')
        for m in self.meas:
            if m.tipo != 'i':
                print(f'{m.tipo}-{str(m.line.nodes[0].ref)+"-"+str(m.line.nodes[1].ref) if hasattr(m,"line") else m.node.ref}: \t {m.value:.3f}  \t {m.h():.3f}')
            else:
                print(f'{m.tipo}-{str(m.line.nodes[0].ref)+"-"+str(m.line.nodes[1].ref) if hasattr(m,"line") else m.node.ref}: \t {np.sqrt(m.value):.3f}  \t {np.sqrt(m.h()):.3f}')
        print('')
         
class node:
    def __init__(self, ref):
        self.ref = ref   
        self.lines = list()
        self.neigh = list()
        self.meas = list()
        self.V = 1
        self.theta = 0
        
    def check(self, currents = 0):
        Ilines = 0
        for line in self.lines:
            if line.nodes[0] == self:
                Ilines += line.I
            else:
                Ilines -= line.I
        Iloads = - complex(self.P, - self.Q)/np.conjugate(self.U) + currents
        return Ilines + Iloads 
    
class line:
    def __init__(self, ref, From, To, R, X, nodes_list):
        self.ref = ref     
        self.Z = complex(R, X)  
        self.G, self.B = -np.real(1/self.Z), -np.imag(1/self.Z)
        self.Y = 1/self.Z
        self.nodes = [next((item for item in nodes_list if item.ref == From), None), 
                      next((item for item in nodes_list if item.ref == To), None)]   
        self.nodes[0].lines.append(self)
        self.nodes[1].lines.append(self)
        self.nodes[0].neigh.append(self.nodes[1])
        self.nodes[1].neigh.append(self.nodes[0])
        self.meas = list()
        
    def check(self):
        res = self.Z*self.I - (self.nodes[0].U - self.nodes[1].U)
        res = self.I - self.Y*(self.nodes[0].U - self.nodes[1].U) 
        return res
  
class measurement:
    def __init__(self, ref, node_id, line_id, tipo, value, std, nodes, lines, n):
        self.ref = ref
        self.tipo = tipo
        self.value = value
        self.std = std
        self.n = n
        if node_id:
            for n in nodes:
                if n.ref == node_id:
                    self.node = n
                    n.meas.append(self)
                    break
        if line_id:
            for l in lines:
                if l.ref == line_id:
                    self.line = l
                    l.meas.append(self)
                    break        
        if hasattr(self, 'node'):
            aux = [[self.node.pointer[0]], [self.node.pointer[1]], 
                  [n.pointer[0] for n in self.node.neigh], 
                  [n.pointer[1] for n in self.node.neigh]]
            self.pointer = [it for item in aux for it in item]
        if hasattr(self, 'line'):
            aux = [[self.line.nodes[0].pointer[0]], [self.line.nodes[0].pointer[1]], 
                  [self.line.nodes[1].pointer[0]], [self.line.nodes[1].pointer[1]]]
            self.pointer = [it for item in aux for it in item]
        
    def h(self):
        if self.tipo == 'p' and hasattr(self, 'node'):
            h = self.Pi()
        if self.tipo == 'q' and hasattr(self, 'node'):
            h = self.Qi()
        if self.tipo == 'i' and hasattr(self, 'node'):
            h = self.Ii()
        if self.tipo == 'p' and hasattr(self, 'line'):
            h = self.Pij()
        if self.tipo == 'q' and hasattr(self, 'line'):
            h = self.Qij()
        if self.tipo == 'i' and hasattr(self, 'line'):
            h = self.Iij()
        if self.tipo == 'v':
            h = self.V()
        return h
        
    def compute_h(self):         
        if self.tipo == 'v' and hasattr(self, 'node'):
            self.dh = np.zeros(self.n)
            self.dh[self.pointer[1]] = 1
        else:
            if self.tipo == 'p' and hasattr(self, 'node'):
                aux = [it for item in [self.Pi_thetai(), self.Pi_Vi(), self.Pi_thetaj(), self.Pi_Vj()] for it in item]        
            if self.tipo == 'q' and hasattr(self, 'node'):
                aux = [it for item in [self.Qi_thetai(), self.Qi_Vi(), self.Qi_thetaj(), self.Qi_Vj()] for it in item]    
            if self.tipo == 'i' and hasattr(self, 'node'):
                aux = [it for item in [self.Ii_thetai(), self.Ii_Vi(), self.Ii_thetaj(), self.Ii_Vj()] for it in item]
            if self.tipo == 'p' and hasattr(self, 'line'):
                aux = [it for item in [self.Pij_thetai(), self.Pij_Vi(), self.Pij_thetaj(), self.Pij_Vj()] for it in item]   
            if self.tipo == 'q' and hasattr(self, 'line'):
                aux = [it for item in [self.Qij_thetai(), self.Qij_Vi(), self.Qij_thetaj(), self.Qij_Vj()] for it in item]  
            if self.tipo == 'i' and hasattr(self, 'line'):
                aux = [it for item in [self.Iij_thetai(), self.Iij_Vi(), self.Iij_thetaj(), self.Iij_Vj()] for it in item]  
                
            self.dh = np.zeros(self.n)
            for item in zip(self.pointer, aux):
                if item[0] != None:
                    self.dh[item[0]] = item[1]
                
    # Derivadas parciales
        
    def Pij_Vi(self, line = None):
        Pij_Vi = list()
        if line == None:
            line = self.line
        node1 = line.nodes[0]
        node2 = line.nodes[1]
        Pij_Vi.append(node2.V*(line.G*np.cos(node1.theta - node2.theta) + line.B*np.sin(node1.theta - node2.theta)) - 2*line.G*node1.V)
        return Pij_Vi
     
    def Pij_Vj(self, line = None):
        Pij_Vj = list()
        if line == None:
            line = self.line
        node1 = line.nodes[0]
        node2 = line.nodes[1]
        Pij_Vj.append(node1.V*(line.G*np.cos(node1.theta - node2.theta) + line.B*np.sin(node1.theta - node2.theta)))
        return Pij_Vj
         
    def Qij_Vi(self, line = None):
        Qij_Vi = list()
        if line == None:
            line = self.line
        node1 = line.nodes[0]
        node2 = line.nodes[1]
        Qij_Vi.append(node2.V*(line.G*np.sin(node1.theta - node2.theta) - line.B*np.cos(node1.theta - node2.theta)) + 2*line.B*node1.V)
        return Qij_Vi
     
    def Qij_Vj(self, line = None):
        Qij_Vj = list()
        if line == None:
            line = self.line
        node1 = line.nodes[0]
        node2 = line.nodes[1]
        Qij_Vj.append(node1.V*(line.G*np.sin(node1.theta - node2.theta) - line.B*np.cos(node1.theta - node2.theta)))
        return Qij_Vj
         
    def Pij_thetai(self, line = None):
        Pij_thetai = list()
        if line == None:
            line = self.line
        node1 = line.nodes[0]
        node2 = line.nodes[1]
        Pij_thetai.append(node1.V*node2.V*(-line.G*np.sin(node1.theta - node2.theta) + line.B*np.cos(node1.theta - node2.theta)))
        return Pij_thetai
     
    def Pij_thetaj(self, line = None):
        Pij_thetaj = list()
        if line == None:
            line = self.line
        node1 = line.nodes[0]
        node2 = line.nodes[1]
        Pij_thetaj.append(node1.V*node2.V*(line.G*np.sin(node1.theta - node2.theta) - line.B*np.cos(node1.theta - node2.theta)))
        return Pij_thetaj
         
    def Qij_thetai(self, line = None):
        Qij_thetai = list()
        if line == None:
            line = self.line
        node1 = line.nodes[0]
        node2 = line.nodes[1]
        Qij_thetai.append(node1.V*node2.V*(line.G*np.cos(node1.theta - node2.theta) + line.B*np.sin(node1.theta - node2.theta)))
        return Qij_thetai
     
    def Qij_thetaj(self, line = None):
        Qij_thetaj = list()
        if line == None:
            line = self.line
        node1 = line.nodes[0]
        node2 = line.nodes[1]
        Qij_thetaj.append(-node1.V*node2.V*(line.G*np.cos(node1.theta - node2.theta) + line.B*np.sin(node1.theta - node2.theta)))
        return Qij_thetaj
    
    def Iij_Vi(self, line = None):
        if line == None:
            line = self.line
        Pij = self.Pij(line)
        Qij = self.Pij(line)
        Vi = line.nodes[0].V
        dPij = self.Pij_Vi(line)[0]
        dQij = self.Qij_Vi(line)[0]
        return [(2/Vi**2)*(Pij*dPij + Qij*dQij) - (2/Vi)*self.value]
        
    def Iij_thetai(self, line = None):
        if line == None:
            line = self.line
        Pij = self.Pij(line)
        Qij = self.Pij(line)
        Vi = line.nodes[0].V
        dPij = self.Pij_thetai(line)[0]
        dQij = self.Qij_thetai(line)[0]
        return [(2/Vi)*(Pij*dPij + Qij*dQij)]
    
    def Iij_Vj(self, line = None):
        if line == None:
            line = self.line
        Pij = self.Pij(line)
        Qij = self.Pij(line)
        Vi = line.nodes[0].V
        dPij = self.Pij_Vj(line)[0]
        dQij = self.Qij_Vj(line)[0]
        return [(2/Vi**2)*(Pij*dPij + Qij*dQij)]
        
    def Iij_thetaj(self, line = None):        
        if line == None:
            line = self.line
        Pij = self.Pij(line)
        Qij = self.Pij(line)
        Vi = line.nodes[0].V
        dPij = self.Pij_thetaj(line)[0]
        dQij = self.Qij_thetaj(line)[0]
        return [(2/Vi)*(Pij*dPij + Qij*dQij)]
    
    def Pi_Vi(self):
        Pi_Vi = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Pi_Vi += self.Pij_Vi(line)
        return [np.sum(Pi_Vi)]
     
    def Pi_Vj(self):
        Pi_Vj = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Pi_Vj += self.Pij_Vj(line)
        return Pi_Vj
         
    def Qi_Vi(self):
        Qi_Vi = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Qi_Vi += self.Qij_Vi(line)
        return [np.sum(Qi_Vi)]
     
    def Qi_Vj(self):
        Qi_Vj = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Qi_Vj += self.Qij_Vj(line)
        return Qi_Vj
         
    def Pi_thetai(self):
        Pi_thetai = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Pi_thetai += self.Pij_thetai(line)
        return [np.sum(Pi_thetai)]
     
    def Pi_thetaj(self):
        Pi_thetaj = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Pi_thetaj += self.Pij_thetaj(line)
        return Pi_thetaj
         
    def Qi_thetai(self):
        Qi_thetai = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Qi_thetai += self.Qij_thetai(line)
        return [np.sum(Qi_thetai)]
     
    def Qi_thetaj(self):
        Qi_thetaj = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Qi_thetaj += self.Qij_thetaj(line)
        return Qi_thetaj
    
    def Ii_Vi(self):        
        Ii_Vi = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Ii_Vi += self.Iij_Vi(line)
        return [np.sum(Ii_Vi)]
        
    def Ii_thetai(self):
        Ii_thetai = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Ii_thetai += self.Iij_thetai(line)
        return [np.sum(Ii_thetai)]
    
    def Ii_Vj(self):   
        Ii_Vj = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Ii_Vj += self.Iij_Vj(line)
        return Ii_Vj
    
    def Ii_thetaj(self):
        Ii_thetaj = list()
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Ii_thetaj += self.Iij_thetaj(line)
        return Ii_thetaj
        
    
    
    # Flujos
        
    def Pi(self):   
        Pi = 0
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Pi += self.node.V*neigh.V*(line.G*np.cos(self.node.theta - neigh.theta) + line.B*np.sin(self.node.theta - neigh.theta))
        return Pi + self.node.V*self.node.V*np.sum([-line.G for line in self.node.lines])
         
    def Qi(self):   
        Qi = 0
        for neigh, line in zip(self.node.neigh, self.node.lines):
            Qi += self.node.V*neigh.V*(line.G*np.sin(self.node.theta - neigh.theta) - line.B*np.cos(self.node.theta - neigh.theta))
        return Qi + self.node.V*self.node.V*np.sum([line.B for line in self.node.lines])
         
    def Ii(self):
        Pi = self.Pi()
        Qi = self.Qi()
        return (Pi**2 + Qi**2) / self.node.V**2
        
    def Pij(self, line = None):         
        if line == None:
            line = self.line
        node1 = line.nodes[0]
        node2 = line.nodes[1]
        Pij = node1.V*node2.V*(line.G*np.cos(node1.theta - node2.theta) + line.B*np.sin(node1.theta - node2.theta)) - line.G*(node1.V**2)
        return Pij
         
    def Qij(self, line = None):        
        if line == None:
            line = self.line
        node1 = line.nodes[0]
        node2 = line.nodes[1]
        Qij = node1.V*node2.V*(line.G*np.sin(node1.theta - node2.theta) - line.B*np.cos(node1.theta - node2.theta)) + line.B*(node1.V**2)
        return Qij        
    
    def Iij(self, line = None):       
        if line == None:
            line = self.line
        Pij = self.Pij()
        Qij = self.Qij()
        return (Pij**2 + Qij**2) / line.nodes[0].V**2
        
    def V(self):
        return self.node.V
        
        
        
        
        
        
        
        
        
    