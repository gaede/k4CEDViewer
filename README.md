# k4CEDViewer

Port of the CEDViewer package from iLCSoft to use Gaudi and EDM4hep (replacing Marlin and LCIO).

### this is WIP ...

## Dependencies

* CED

* PODIO/EDM4hep

* Gaudi

* k4FWCore



### Compilation

Run, from the `k4CEDViewer` directory:

``` bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh
k4_local_repo
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install -G Ninja
ninja install
```

Alternatively you can source the nightlies instead of the releases:

``` bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```

Note that if you source the releases and use the current version of this
repository this is not guaranteed to work as there could be changes since this
repository was built for the release. What you can do in this case is to
checkout a previous tag, for example:

``` bash
git checkout v0.5.0
```

This is because the releases are only built with tagged versions of the
packages. With the nightlies this repository should always work; if it doesn't
please [open an
issue](https://github.com/key4hep/k4CEDViewer/issues/new/choose).

### Execute Examples

Make sure that `k4CEDViewer/install/lib` and
`k4CEDViewer/install/python` are in `LD_LIBRARY_PATH` and `PYTHONPATH`
respectively (`k4_local_repo` should take care of this). If they are not, they
can be added by running:
``` bash
export LD_LIBRARY_PATH=$PWD/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/install/python:$PYTHONPATH
```

and then run the examples like this:


``` bash
glced &
k4run ./k4CEDViewer/options/k4_event_display.py --Input=../ILDConfig/StandardConfig/production/E250-SetA.P4f_ww_h.Gwhizard-2_8_5.eL.pR.I500066.1_REC.edm4hep.root --compactFile=$k4geo_DIR/FCCee/ILD_FCCee/compact/ILD_FCCee_v01/ILD_FCCee_v01.xml --log-level=info
```


### Continuous integration (CI) in forks
If this repository is forked instead of used as a template, CI will not work by
default, because GitHub disables workflows by default in forks of foreign
repositories. Workflows can be enabled by going to the Actions tab and clicking
the green button labeled `I understand my workflows, go ahead and enable them`.
If this is not done, no workflows will run, for example, after making a pull
request. In addition, scheduled workflows are also disabled by default and
require an additional step to be enabled manually from the Actions tab. The
workflows that are disabled will have `disabled` next to them (at the time of
writing, the only one is the `Key4hep build`).


## References:
These could perhaps be usefule for newcomers.
1. [lhcb-98-064 COMP](https://cds.cern.ch/record/691746/files/lhcb-98-064.pdf)
2. [Hello World in the Gaudi Framework](https://lhcb.github.io/DevelopKit/02a-gaudi-helloworld)
