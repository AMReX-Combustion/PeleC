import yt
import argparse


parser = argparse.ArgumentParser(
            description = 'Post Processing script to take field average of a variable.'
            )
parser.add_argument(
            'plotfile',
            type = str,
            help = 'First argument must be the name of the plotfile.'
            )
parser.add_argument(
            'name_var',
            type = str,
            help = 'Name of the variable: x_velocity, MachNumber, Temp, density, ...'
            )
args = parser.parse_args()
# Remove Trailing Slashes
args.plotfile = '/'.join(list(filter(lambda x: len(x), args.plotfile.split('/'))))

#print("{}".format(args))


ds = yt.load(args.plotfile)  # load data

field = args.name_var  # The field to average
weight = "cell_mass"  # The weight for the average

ad = ds.all_data()  # This is a region describing the entire box,
                    # but note it doesn't read anything in yet!

# We now use our 'quantities' call to get the average quantity
average_value = ad.quantities.weighted_average_quantity(field, weight)

print("Average %s (weighted by %s) is %0.3e %s" % (field, weight, average_value, average_value.units))
