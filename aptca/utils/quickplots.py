import matplotlib.pyplot as plt



def plot_density_distribution(Density,distance_axial,density_axial,
                                distance_radial,density_radial,save=False,save_path='.',save_name='density_distr'):
    """
    plot result from AtomCloud.get_atom_density_distribution() function
    """
    fig, ax = plt.subplots(3,1,figsize=[8,8])
    #fig.suptitle("Atom Density Distributions",fontsize=20)
    # plot axial density distribution
    ax[0].plot(distance_axial,density_axial,color='k',linewidth=2)
    ax[0].set_ylabel('Density ($\mathregular{at/nm^3}$)')
    ax[0].set_xlabel('Depth (nm)')
    ax[0].set_title('Axial Density Distribution',fontsize=16)
    # radial
    ax[1].plot(distance_radial,density_radial,color='k',linewidth=2)
    ax[1].set_ylabel('Density ($\mathregular{at/nm^3}$)')
    ax[1].set_xlabel('Radius (nm)')
    ax[1].set_title('Radial Density Distribution',fontsize=16)
    # histogram
    ax[2].hist(Density,bins=100)
    ax[2].set_ylabel('Counts')
    ax[2].set_xlabel('Density ($\mathregular{at/nm^3}$)')
    ax[2].set_title('Density Histogram',fontsize=16)
    fig.tight_layout()
    if save:
        fig.savefig(save_path+'/'+save_name+'.jpg',dpi=300,bbox_inches='tight')