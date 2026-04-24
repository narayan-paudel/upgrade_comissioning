class xDOMSensorMeasurement:
    def __init__(self,df,string,port):
        self.df = df
        self.string = string
        self.port = port
        self.get_mean_B()
        self.get_mean_g()

    def get_mean_B(self):
        '''
        returns field in muT and degrees
        '''
        bx,by,bz = self.df[["bx","by","bz"]].values.T
        B,B_theta,B_phi = to_spherical_list(bx,by,bz)
        if np.rad2deg(np.mean(B_phi))<0:
            print(f"negative  phi detected {np.rad2deg(np.mean(B_phi))}")
        for ielt,iname in zip([bx,by,bz,B],["bx","by","bz","B"]):
            setattr(self, f"{iname}_mean", np.mean(ielt)*10**6)
            setattr(self, f"{iname}_std", np.std(ielt)*10**6)
        for ielt,iname in zip([B_theta,B_phi],["B_theta","B_phi"]):
            setattr(self, f"{iname}_mean", np.rad2deg(np.mean(ielt)))
            setattr(self, f"{iname}_std", np.rad2deg(np.std(ielt)))


    def get_mean_g(self):
        '''
        returns field in ms^-2 and degrees
        '''
        gx,gy,gz = self.df[["gx","gy","gz"]].values.T
        g,g_theta,g_phi = to_spherical_list(gx,gy,gz)
        if np.rad2deg(np.mean(g_phi))<0:
            print(f"negative  phi detected {np.rad2deg(np.mean(g_phi))}")
        for ielt,iname in zip([gx,gy,gz,g],["gx","gy","gz","g"]):
            setattr(self, f"{iname}_mean", np.mean(ielt))
            setattr(self, f"{iname}_std", np.std(ielt))
        for ielt,iname in zip([g_theta,g_phi],["g_theta","g_phi"]):
            setattr(self, f"{iname}_mean", np.rad2deg(np.mean(ielt)))
            setattr(self, f"{iname}_std", np.rad2deg(np.std(ielt)))