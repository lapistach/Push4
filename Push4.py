import numpy as np
from matplotlib import cm
from matplotlib.collections import EllipseCollection
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pickle
            
class CollisionEvent:
 
    def __init__(self, Type = 'wall or other', dt = np.inf, mono_1 = 0, mono_2 = 0, w_dir = 1):
      
        self.Type = Type
        self.dt = dt
        self.mono_1 = mono_1
        self.mono_2 = mono_2  # only importent for interparticle collisions
        self.w_dir = w_dir # only important for wall collisions
        
        
    def __str__(self):
        if self.Type == 'wall':
            return "Event type: {:s}, dt: {:.8f}, p1 = {:d}, dim = {:d}".format(self.Type, self.dt, self.mono_1, self.w_dir)
        else:
            return "Event type: {:s}, dt: {:.8f}, p1 = {:d}, p2 = {:d}".format(self.Type, self.dt, self.mono_1, self.mono_2)

class Monomers:

    def __init__(self, NumberOfMonomers = 4, L_xMin = 0, L_xMax = 1, L_yMin = 0, L_yMax = 1, NumberMono_per_kind = np.array([4]), Radiai_per_kind = 0.5*np.ones(1), Densities_per_kind = np.ones(1), k_BT = 1, FilePath = './Configuration.p'):
        try:
            self.__dict__ = pickle.load( open( FilePath, "rb" ) )
            print("IMPORTANT! System is initialized from file %s, i.e. other input parameters of __init__ are ignored!" % FilePath)
        except:
            assert ( NumberOfMonomers > 0 )
            assert ( (L_xMin < L_xMax) and (L_yMin < L_yMax) )
            self.NM = NumberOfMonomers
            self.DIM = 2 #dimension of system
            self.BoxLimMin = np.array([ L_xMin, L_yMin])
            self.BoxLimMax = np.array([ L_xMax, L_yMax])
            self.mass = -1*np.ones( self.NM ) # Masses, negative mass means not initialized
            self.rad = -1*np.ones( self.NM ) # Radiai, negative radiai means not initialized
            self.pos = np.empty( (self.NM, self.DIM) ) # Positions, not initalized but desired shape
            self.vel = np.empty( (self.NM, self.DIM) ) # Velocities, not initalized but desired shape
            self.mono_pairs = np.array( [ (k,l) for k in range(self.NM) for l in range( k+1,self.NM ) ] )
            self.next_wall_coll = CollisionEvent( 'wall', np.inf, 0, 0, 0)
            self.next_mono_coll = CollisionEvent( 'mono', np.inf, 0, 0, 0)
            self.p= 'Planetesimal is forming!'
            self.assignRadiaiMassesVelocities(NumberMono_per_kind, Radiai_per_kind, Densities_per_kind, k_BT )
            self.assignRandomMonoPos( )
    
    def save_configuration(self, FilePath = 'MonomerConfiguration.p'):
        '''Saves configuration. Callable at any time during simulation.'''
        #print( self.__dict__ )
    
    def assignRadiaiMassesVelocities(self, NumberMono_per_kind = np.array([4]), Radiai_per_kind = 0.5*np.ones(1), Densities_per_kind = np.ones(1), k_BT = 1 ):
        self.vel[self.NM-1] = [0,0]
        assert( sum(NumberMono_per_kind) == self.NM )
        assert( isinstance(Radiai_per_kind,np.ndarray) and (Radiai_per_kind.ndim == 1) )
        assert( (Radiai_per_kind.shape == NumberMono_per_kind.shape) and (Radiai_per_kind.shape == Densities_per_kind.shape))
  
        total=0
        for i, number in enumerate(NumberMono_per_kind):    
          self.rad [total:total+number]=Radiai_per_kind[i]
          self.mass[total:total+number]=Densities_per_kind[i]*(np.pi*(Radiai_per_kind[i]))
          total=number

        assert( k_BT > 0 )
        self.vel = np.random.rand(self.NM,self.DIM)
        c = 0
        for i in range(self.NM-1):
            c += (self.mass[i]/2)*(self.vel[i,0]**2+self.vel[i,1]**2)
        g = (self.NM-1)*self.DIM*k_BT/c
        self.vel *= np.ones(np.shape(self.vel))*np.sqrt(g)
        self.vel[self.NM-1]=[0,0]
        
    
    def assignRandomMonoPos(self, start_index = 0 ):

        assert ( min(self.rad) > 0 )#otherwise not initialized
        mono_new, infiniteLoopTest = start_index, 0
        BoxLength = self.BoxLimMax - self.BoxLimMin
        while mono_new < self.NM and infiniteLoopTest < 10**4:
            infiniteLoopTest += 1
            self.pos[self.NM - 1, :] = [(L_xMax-L_xMin)/2,(L_yMax-L_yMin)/2]
           
    
    def __str__(self, index = 'all'):
        if index == 'all':
            return "\nMonomers with:\nposition = " + str(self.pos) + "\nvelocity = " + str(self.vel) + "\nradius = " + str(self.rad) + "\nmass = " + str(self.mass)
        else:
            return "\nMonomer at index = " + str(index) + " with:\nposition = " + str(self.pos[index]) + "\nvelocity = " + str(self.vel[index]) + "\nradius = " + str(self.rad[index]) + "\nmass = " + str(self.mass[index])
        
    def Wall_time(self):

        coll_condition = np.where( self.vel > 0, self.BoxLimMax-self.rad[:,np.newaxis], self.BoxLimMin+self.rad[:,np.newaxis])
        dt_List = ( coll_condition - self.pos) / self.vel
        MinTimeIndex = np.argmin(  dt_List )
        collision_disk = MinTimeIndex // 2
        wall_direction = MinTimeIndex % 2
        
        self.next_wall_coll.dt = dt_List[collision_disk][wall_direction]
        self.next_wall_coll.mono_1 = collision_disk
        self.next_wall_coll.w_dir = wall_direction
        
        print(self.next_wall_coll)
    
        
    def Mono_pair_time(self):
       

        mono_i = self.mono_pairs[:,0] # List of collision partner 1
        mono_j = self.mono_pairs[:,1] # List of collision partner 2
        
        dx_dy_ofPairs = self.pos[mono_i] - self.pos[mono_j]

        vx_vy_ofPairs= self.vel[mono_i]-self.vel[mono_j]

        print(vx_vy_ofPairs)

        a=vx_vy_ofPairs[:,0]**2+vx_vy_ofPairs[:,1]**2
        b=2*(vx_vy_ofPairs[:,0]*dx_dy_ofPairs[0,0]+vx_vy_ofPairs[:,1]*dx_dy_ofPairs[0,1])
        c=(dx_dy_ofPairs[0,0])**2+ (dx_dy_ofPairs[0,1])**2 - (self.rad[mono_i]+self.rad[mono_i])**2
        delta=b**2-4*a*c

        deltat= np.where(((delta>0)&(b<0)), (1/(2*a))*(-b-np.sqrt(delta)), np.inf)
        
        min_index=np.argmin(deltat)
        
        collision_disk_1, collision_disk_2 = self.mono_pairs[min_index]
        
        
        self.next_mono_coll.dt = deltat[min_index]
        self.next_mono_coll.mono_1 = collision_disk_1
        self.next_mono_coll.mono_2 = collision_disk_2

        print(self.next_mono_coll)

    def compute_next_event(self):
        
        self.Wall_time()
        self.Mono_pair_time()
        if self.next_wall_coll.dt < self.next_mono_coll.dt :
            return self.next_wall_coll
        else :
            return self.next_mono_coll

            
    def compute_new_velocities(self, next_event):
        
        if next_event.Type=='wall':

            self.vel[next_event.mono_1,next_event.w_dir]= self.vel[next_event.mono_1,next_event.w_dir]*-1

        else:
            mono_1 = next_event.mono_1
            mono_2 = next_event.mono_2
            diff = self.pos[mono_2] - self.pos[mono_1]
            diff=diff/np.linalg.norm(diff)
            m1=self.mass[mono_1]
            m2=self.mass[mono_2]
            Dv=self.vel[mono_1]-self.vel[mono_2]
            self.vel[mono_1]=self.vel[mono_1] - (2*m2)/(m1+m2)*np.inner(diff,Dv)*diff
            self.vel[mono_2]=self.vel[mono_2] + (2*m1)/(m1+m2)*np.inner(diff,Dv)*diff
           
            
            NumberOfMobileMonomers= self.NM
            for i in range(self.NM):
                if self.vel[i][0]==0:
                   NumberOfMobileMonomers -= 1

            if NumberOfMobileMonomers == 1:

                if self.vel[next_event.mono_1][0]==0:

                   Vd=self.vel[next_event.mono_2]
                   rapport=diff[1]/diff[0]
                   new_x=-2*(rapport*Vd[1]+Vd[0])/(1+rapport**2)
                   new_y=rapport*new_x
                   new=[new_x,new_y]
                   self.vel[next_event.mono_2]+=new
              
                else :

                    Vd=self.vel[next_event.mono_1]
                    rapport=diff[1]/diff[0]
                    new_x=-2*(rapport*Vd[1]+Vd[0])/(1+rapport**2)
                    new_y=rapport*new_x
                    new=[new_x,new_y]
                    self.vel[next_event.mono_1]+=new


            
                   
    def snapshot(self, FileName = './snapshot.png', Title = '$t = $?'):

        fig, ax = plt.subplots( dpi=300 )
        L_xMin, L_xMax = self.BoxLimMin[0], self.BoxLimMax[0]
        L_yMin, L_yMax = self.BoxLimMin[1], self.BoxLimMax[1]
        BorderGap = 0.1*(L_xMax - L_xMin)
        ax.set_xlim(L_xMin-BorderGap, L_xMax+BorderGap)
        ax.set_ylim(L_yMin-BorderGap, L_yMax+BorderGap)

        #--->plot hard walls (rectangle)
        rect = mpatches.Rectangle((L_xMin,L_yMin), L_xMax-L_xMin, L_yMax-L_yMin, linestyle='dashed', ec='gray', fc='None')
        ax.add_patch(rect)
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ position')
        ax.set_ylabel('$y$ position')
        
        #--->plot monomer positions as circles
        MonomerColors = np.linspace( 0.2, 0.95, self.NM)
        Width, Hight, Angle = 2*self.rad, 2*self.rad, np.zeros( self.NM )
        collection = EllipseCollection( Width, Hight, Angle, units='x', offsets=self.pos,
                       transOffset=ax.transData, cmap='nipy_spectral', edgecolor = 'k')
        collection.set_array(MonomerColors)
        collection.set_clim(0, 1) # <--- we set the limit for the color code
        ax.add_collection(collection)

        #--->plot velocities as arrows
        ax.quiver( self.pos[:,0], self.pos[:,1], self.vel[:,0], self.vel[:,1] , units = 'dots', scale_units = 'dots')
        
        plt.title(Title)
        plt.savefig(FileName)
        plt.close()


import os
from matplotlib import cm
from matplotlib.collections import EllipseCollection
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

np.random.seed(999)
Snapshot_output_dir = './SnapshotsMonomers'
if not os.path.exists(Snapshot_output_dir): os.makedirs(Snapshot_output_dir)

Conf_output_dir = './ConfsMonomers'
Path_ToConfiguration = Conf_output_dir+'/FinalMonomerConf.p'
if False: #os.path.isfile( Path_ToConfiguration ):
    '''Initialize system from existing file'''
    mols = Monomers( FilePath = Path_ToConfiguration )
else:
    '''Initialize system with following parameters'''
    # create directory if it does not exist
    if not os.path.exists(Conf_output_dir): os.makedirs(Conf_output_dir)
    #define parameters
    NumberOfMonomers = 13
    L_xMin, L_xMax = 0, 100
    L_yMin, L_yMax = 0, 50
    NumberMono_per_kind = np.array([NumberOfMonomers-1,1]) 
    Radiai_per_kind = np.array([ 1.5,1.5]) 
    Densities_per_kind = np.append(np.ones(len(NumberMono_per_kind)-1),np.inf)
    k_BT = 1000
    # call constructor, which should initialize the configuration
    mols = Monomers(NumberOfMonomers, L_xMin, L_xMax, L_yMin, L_yMax, NumberMono_per_kind, Radiai_per_kind, Densities_per_kind, k_BT )
    
mols.snapshot( FileName = Snapshot_output_dir+'/InitialConf.png', Title = '$t = 0$')
#we could initialize next_event, but it's not necessary
#next_event = pc.CollisionEvent( Type = 'wall or other, to be determined', dt = 0, mono_1 = 0, mono_2 = 0, w_dir = 0)

t = 0.0
dt = 0.02
NumberOfFrames = 1000
next_event = mols.compute_next_event()
#next_time_vel_change=t+next_event.dt
def MolecularDynamicsLoop( frame ):

    global t, mols, next_event,next_time_vel_change
    t+=dt
    while t >= next_event.dt:
        mols.compute_new_velocities(next_event)
        next_event=mols.compute_next_event()
        next_event.dt += t
    mols.pos += mols.vel * dt   
            
   # timeremaining=future_time_next_frame-t
   # mols.pos += timeremaining*mols.vel
 #   t += timeremaining
 #   next_event.dt -= timeremaining
    
    plt.title( '$t = %.4f$, remaining frames = %d' % (t, NumberOfFrames-(frame+1)) )
    collection.set_offsets( mols.pos )
    return collection


fig, ax = plt.subplots()
L_xMin, L_yMin = mols.BoxLimMin #not defined if initalized by file
L_xMax, L_yMax = mols.BoxLimMax #not defined if initalized by file
BorderGap = 0.1*(L_xMax - L_xMin)
ax.set_xlim(L_xMin-BorderGap, L_xMax+BorderGap)
ax.set_ylim(L_yMin-BorderGap, L_yMax+BorderGap)
ax.set_aspect('equal')

# confining hard walls plotted as dashed lines
rect = mpatches.Rectangle((L_xMin,L_yMin), L_xMax-L_xMin, L_yMax-L_yMin, linestyle='dashed', ec='gray', fc='None')
ax.add_patch(rect)


# plotting all monomers as solid circles of individual color
MonomerColors = np.linspace(0.2,0.95,mols.NM)
Width, Hight, Angle = 2*mols.rad, 2*mols.rad, np.zeros(mols.NM)
collection = EllipseCollection(Width, Hight, Angle, units='x', offsets=mols.pos,
                       transOffset=ax.transData, cmap='nipy_spectral', edgecolor = 'k')
collection.set_array(MonomerColors)
collection.set_clim(0, 1) # <--- we set the limit for the color code
ax.add_collection(collection)

'''Create the animation, i.e. looping NumberOfFrames over the update function'''
Delay_in_ms = 33.3 #delay between images/frames for plt.show()
ani = FuncAnimation(fig, MolecularDynamicsLoop, frames=NumberOfFrames, interval=Delay_in_ms, blit=False, repeat=False)
plt.show()

'''Save the final configuration and make a snapshot.'''
#write the function to save the final configuration
mols.save_configuration(Path_ToConfiguration)
mols.snapshot( FileName = Snapshot_output_dir + '/FinalConf.png', Title = '$t = %.4f$' % t)
