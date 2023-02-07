# Handling optical link data

`tintervals.rocitlinks` provides a structure to handle frequency data in the format established in the [EMPIR project ROCIT](http://empir.npl.co.uk/rocit/).

The package can handle data in fractional frequency or transfer beat notation based on the math in
[Lodewyck et al.,  "Universal formalism for data sharing and processing in clock comparison networks", Phys. Rev. Research  2, 043269 (2020).](https://link.aps.org/doi/10.1103/PhysRevResearch.2.043269)


## Basic usage

    import tintervals.rocitlinks as rl

    # load link data
    link1 = rl.load_link_from_dir('./Data/Link1', meta='metadata.yml')
    link2 = rl.load_link_from_dir('./Data/Link2', meta='metadata.yml')

    # chaining link
    reslink, masks = rl.chain(link1, link2)

    # averaging links
    days, daylink, dcount = rl.link_average(reslink, 'day')

    # saving resulting link
    rl.save_link_to_dir('./Output', daylink)


## Module documentation
```{eval-rst}
.. automodule:: tintervals.rocitlinks
  :members:
  :autosummary:
```