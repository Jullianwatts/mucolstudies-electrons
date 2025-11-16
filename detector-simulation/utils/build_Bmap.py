import ROOT
import math
from array import array


# Function to calculate the magnetic field components Bx, By for a given (x, y) point
def calculate_magnetic_field(x, y, magnitude):
    # Calculate the angle theta of the point (x, y) with respect to the positive x-axis
    theta = math.atan2(y, x)

    # Calculate the magnetic field components Bx, By
    Bx = magnitude * math.cos(theta)
    By = magnitude * math.sin(theta)

    return Bx, By


def create_and_fill_tree(Bfield, rmin, rmax, halflength, z_endcap_min, num_pointsz, num_pointsxy):
    # Create a TTree
    tree = ROOT.TTree("my_magnetic_field_map",
                      "MuColl_10TeV_v0A toroidal magnetic field")

    # Define variables for the branches
    x_mm = array('f', [0])
    y_mm = array('f', [0])
    z_mm = array('f', [0])
    Bx = array('f', [0])
    By = array('f', [0])
    Bz = array('f', [0])

    # Add branches to the TTree
    tree.Branch('x_mm', x_mm, 'var/F')
    tree.Branch('y_mm', y_mm, 'var/F')
    tree.Branch('z_mm', z_mm, 'var/F')
    tree.Branch('Bx', Bx, 'var/F')
    tree.Branch('By', By, 'var/F')
    tree.Branch('Bz', Bz, 'var/F')

    for iz in range(num_pointsz+1):
        z = -halflength + (iz * 2. * halflength) / num_pointsz

        print("z = ", z)

        for ix in range(num_pointsxy+1):
            x = -rmax + (ix * 2. * rmax) / num_pointsxy

            for iy in range(num_pointsxy+1):

                y = -rmax + (iy * 2. * rmax) / num_pointsxy
                radius = math.sqrt(x*x+y*y)

                x_mm[0] = x
                y_mm[0] = y
                z_mm[0] = z
                Bx_val, By_val = calculate_magnetic_field(x, y, Bfield)
                if radius < rmin and math.fabs(z) < z_endcap_min:
                    Bx_val = 0.
                    By_val = 0.
                Bx[0] = Bx_val
                By[0] = By_val
                Bz[0] = 0.0

                tree.Fill()

    return tree


if __name__ == "__main__":
    # Set the parameters for the cylindrical surface
    rmin = 4090  # min radius [mm]
    rmax = 6100  # max radius [mm]
    halflength = 6398  # Chosen halflength [mm]
    z_endcap_min = 4570  # End of HCal [mm]
    num_pointsz = 28  # Number of points to create on z
    num_pointsxy = 500  # Number of points to create (same for x,y)
    Bfield = 5  # [T]

    # Create and fill the TTree
    tree = create_and_fill_tree(
        Bfield, rmin, rmax, halflength, z_endcap_min, num_pointsz, num_pointsxy)

    # Create a ROOT file to store the TTree
    output_file = ROOT.TFile("cylindrical_surface_data.root", "RECREATE")
    output_file.cd()

    # Write the TTree to the ROOT file
    tree.Write()

    # Close the ROOT file
    output_file.Close()

    print("TTree has been created and filled with data and saved to 'cylindrical_surface_data.root'.")
