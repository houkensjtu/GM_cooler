import taichi as ti

ti.init()

@ti.data_oriented
class max_pressure:
    def __init__(self, dt, t_end, init_temp, copper_mass, heat_load):
        self.dt = dt
        self.t_end = t_end
        self.nstep = int(t_end / dt)
        
        self.init_temp = init_temp
        self.copper_mass = copper_mass
        self.heat_load = heat_load
        
        self.T = ti.field(ti.f32, shape = self.nstep)

        self.cold_T = ti.field(ti.f32, shape = self.nstep)
        self.rege_T = ti.field(ti.f32, shape = self.nstep)
        self.warm_T = ti.field(ti.f32, shape = self.nstep)

        self.cold_V = 0.00045216
        self.rege_V = 0.00025622
        self.warm_V = 9.129e-5

        self.hous_V = 0.00043325

        self.mdot = 1.4197e-6

        self.P_cold = ti.field(ti.f32, shape = (self.nstep,2))
        self.P_warm = ti.field(ti.f32, shape = (self.nstep,2))        

        self.m_cold = ti.field(ti.f32, shape = self.nstep)
        self.m_warm = ti.field(ti.f32, shape = self.nstep)                        


    def log10(self,b):
        return ti.log(b)/ti.log(10)
    

    def copper_specheat(self, t):
        a = -1.91844
        b = -0.15973
        c = 8.61013
        d = -18.996
        e = 21.9661
        f = -12.7328
        g = 3.54322
        h = -0.3797
        
        logy = a + b*self.log10(t) + c*(self.log10(t))**2 + d*(self.log10(t))**3 +\
        e*(self.log10(t))**4 + f*(self.log10(t))**5 + g*(self.log10(t))**6 + h*(self.log10(t))**7
        
        return 10**logy
    

    def solve_T(self):
        for i in range(self.nstep):
            self.T[i] = self.init_temp

        for i in range(1,self.nstep):
            heat_load = (300 - self.T[i-1]) / (300 - self.init_temp) * self.heat_load
            self.T[i] = self.T[i-1] + heat_load * self.dt / (self.copper_mass * self.copper_specheat(self.T[i-1]))

            
    def init_P(self):
        self.m_cold[0] = 0.00846468
        self.P_cold[0, 0] = 2.3

        self.m_warm[0] = 0.7 * self.hous_V * 10**6  / 2077.1 / 300
        self.P_warm[0,0] = 0.7
        
        self.cold_T[0] = self.init_temp
        self.warm_T[0] = 300.0
        self.rege_T[0] = 0.5 * (self.warm_T[0] + self.cold_T[0])        
        for i in range(1,self.nstep):
            self.P_cold[i,0] = 2.3
            self.P_warm[i,0] = 0.7
            self.warm_T[i] = 300.0
            self.cold_T[i] = self.T[i]            
            self.rege_T[i] = 0.5 * (self.warm_T[i] + self.cold_T[i])


    def solve_P(self):
        self.init_P()
        for j in range(1):
            for i in range(1, self.nstep):
                mdot = 4.9743e-6 * ti.exp(4.5424 * (self.P_cold[i-1, 0] - self.P_warm[i-1,0])) * 0.001 * 1.1179
                self.m_cold[i] = self.m_cold[i-1] - mdot * self.dt
                self.m_warm[i] = self.m_warm[i-1] + mdot * self.dt                
            
                self.warm_T[i] = 300.0
                self.cold_T[i] = self.T[i]
                self.rege_T[i] = 0.5 * (self.warm_T[i] + self.cold_T[i])

                RT1 = self.warm_V / (2077.1 * self.warm_T[i])
                RT2 = self.cold_V / (2077.1 * self.cold_T[i])
                RT3 = self.rege_V / (2077.1 * self.rege_T[i])

                RT4 = self.hous_V / (2077.1 * 300)

                self.P_cold[i,0] = self.m_cold[i] / (RT1 + RT2 + RT3) / 10**6
                self.P_warm[i,0] = self.m_warm[i] / RT4 / 10**6
                if self.P_warm[i,0] > 1.93:
                    self.P_warm[i,0] = 1.93
                # self.P_warm[i, 0] = 0.7

                
    def write_data(self, timestep):
        import csv
        header = ['t', 'cold_T', 'regen_T', 'warm_T', 'm', 'P_cold_old', 'P_cold', 'P_warm']
        with open('data.csv', 'w') as f:
            f_csv = csv.writer(f)
            f_csv.writerow(header)
            for i in range(0, self.nstep, int(timestep / self.dt)):
                f_csv.writerow([self.dt*i, self.cold_T[i], self.rege_T[i], self.warm_T[i], self.m_cold[i], self.P_cold[i,0], self.P_cold[i,1], self.P_warm[i,0]])

                
    def solve(self):
        self.solve_T()
        self.solve_P()
        self.write_data(1)
        

def main():
    ch160 = max_pressure(1, 10800, 80, 6.78, 600)
    ch160.solve()

    
if __name__ == "__main__":
    main()
