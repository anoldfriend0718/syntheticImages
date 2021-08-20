import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import subprocess
import os
import shutil
import json
import mpslib as mps


class SedimentBearingRockGenerator(object):
    def __init__(self) -> None:
        super().__init__()
        self._label_pore = 1
        self._label_rock = 2
        self._label_sediment = 3

    def simulate(self, parameters):
        parameters["work_dir"] = os.path.abspath(parameters["work_dir"])
        parameters["rock_simulation"]["simulation_rock_file"] = os.path.abspath(
            parameters["rock_simulation"]["simulation_rock_file"])
        ti_nx = parameters["rock_simulation"]["training_rock_nx"]
        ti_ny = parameters["rock_simulation"]["training_rock_ny"]
        training_rock_file = parameters["rock_simulation"]["training_rock_file"]
        with open(training_rock_file, "r") as fp:
            raw_data = np.fromfile(fp, dtype=np.uint8)
        ti = raw_data.reshape(ti_nx, ti_ny)
        print("training rock shape: ", ti.shape)
        ti_filename = "ti.dat"
        mps.eas.write_mat(ti, ti_filename)

        method = parameters["rock_simulation"]["simulation_method"]
        mps_imp = mps.mpslib(method=method, verbose_level=-1)

        mps_imp.par['n_cond'] = parameters["rock_simulation"]["n_cond"]  # number of the conditional points
        mps_imp.par['template_size'] = np.array(json.loads(parameters["rock_simulation"]["template_size"]))
        mps_imp.par["n_multiple_grids"] = parameters["rock_simulation"]["n_multiple_grids"]
        mps_imp.par['n_real'] = 1  # number of the realization number
        mps_imp.par['n_threads'] = 1  # number of the threads for parallization
        mps_imp.par['do_entropy'] = 1
        si_nx = parameters["rock_simulation"]["simulation_rock_nx"]
        si_ny = parameters["rock_simulation"]["simulation_rock_ny"]
        mps_imp.par['simulation_grid_size'] = np.array([si_nx, si_ny, 1])
        mps_imp.par['ti_fnam'] = ti_filename
        print(mps_imp.par)
        mps_imp.run_parallel()

        simulation_rock = mps_imp.sim[0].reshape(si_nx, si_ny)
        simulation_rock_file = parameters["rock_simulation"]["simulation_rock_file"]
        self._make_parent_dir(simulation_rock_file)
        with open(simulation_rock_file, "w") as fp:
            simulation_rock.astype(np.uint8).tofile(fp)
        ax, fig = self._show_rock(simulation_rock, cmap=colors.ListedColormap(['black', 'gray']))
        img_dir = self._get_img_dir(parameters)
        fig.savefig(f"{img_dir}/simulation_rock.jpg", bbox_inches='tight', dpi=600)

        # clean
        is_clean = parameters["is_clean"]
        if is_clean:
            os.remove(ti_filename)
            os.remove("mps_snesim_tree_000.txt")
            os.remove("ti_thread_000.dat")
            os.remove("ti_thread_000.dat.gslib")
            shutil.rmtree("thread000")

            print("temporatory simulation directories were cleaned")

        return simulation_rock

    def _make_parent_dir(self, file):
        parent = os.path.abspath(os.path.join(file, os.pardir))
        if not os.path.exists(parent):
            os.makedirs(parent, exist_ok=True)

    def deposit(self, simulation_rock, parameters):
        parameters["work_dir"] = os.path.abspath(parameters["work_dir"])
        parameters["sediment_deposition"]["sediment_bearing_rock_file"] = os.path.abspath(
            parameters["sediment_deposition"]["sediment_bearing_rock_file"])
        # transform simulation rock from [0,n] to [0,255]
        max_value = np.max(simulation_rock)
        simulation_rock = 255*(simulation_rock/max_value).astype(np.uint8)

        # random deposit sediment inside the rock
        chunked_subdomain_nx, chunked_subdomain_ny, subdomains_dir, subdomains = self._chunk_rock(
            parameters, simulation_rock)
        self._deposit_sediment_in_subdomains(parameters, chunked_subdomain_nx,
                                             chunked_subdomain_ny, subdomains_dir, subdomains)
        sediment_bearing_subdomains = self._select_random_deposition_subdomains(
            parameters, chunked_subdomain_nx, chunked_subdomain_ny, subdomains_dir)
        sediment_bearing_domain = np.bmat(sediment_bearing_subdomains)
        # output
        self._report_fraction(sediment_bearing_domain)
        sediment_bearing_rock_file = parameters["sediment_deposition"]["sediment_bearing_rock_file"]
        self._make_parent_dir(sediment_bearing_rock_file)
        with open(sediment_bearing_rock_file, "w") as fp:
            sediment_bearing_domain.astype(np.uint8).tofile(fp)
        ax, fig = self._show_rock(sediment_bearing_domain, cmap=colors.ListedColormap(['black', 'gray', 'white']))
        img_dir = self._get_img_dir(parameters)
        fig.savefig(f"{img_dir}/sediment_bearing_rock.jpg", bbox_inches='tight', dpi=600)
        # clean
        is_clean = parameters["is_clean"]
        if is_clean:
            shutil.rmtree(subdomains_dir)
            print("temporatory directories were cleaned")
        return sediment_bearing_domain

    def _get_img_dir(self, parameters):
        img_dir = os.path.join(parameters["work_dir"], "images")
        if not os.path.exists(img_dir):
            os.makedirs(img_dir, exist_ok=True)
        return img_dir

    def write_in_openfoam_format(self, data, parameters):
        parameters["work_dir"] = os.path.abspath(parameters["work_dir"])
        parameters["openfoam_file"]["openfoam_file_output_dir"] = os.path.abspath(
            parameters["openfoam_file"]["openfoam_file_output_dir"])
        parameters["openfoam_file"]["openfoam_template_dir"] = os.path.abspath(
            parameters["openfoam_file"]["openfoam_template_dir"])
        data_2d = data.reshape([data.shape[0], data.shape[1]])
        eps = np.ones([data_2d.shape[0], data_2d.shape[1]])
        eps[:, :] = data_2d[:, :]
        eps[data_2d == self._label_pore] = 1.0
        eps[data_2d == self._label_sediment] = parameters["openfoam_file"]["porosity_at_sediment_region"]
        eps[data_2d == self._label_rock] = parameters["openfoam_file"]["porosity_at_rock_region"]

        sedimemt = np.ones([data_2d.shape[0], data_2d.shape[1]])
        sedimemt[:, :] = data_2d[:, :]
        sedimemt[data_2d == self._label_pore] = 0.0
        sedimemt[data_2d == self._label_sediment] = 1-parameters["openfoam_file"]["porosity_at_sediment_region"]
        sedimemt[data_2d == self._label_rock] = 0.0

        fig, ax = plt.subplots(ncols=2)
        ax[0].imshow(eps, interpolation="none")
        ax[0].set_title("eps")
        ax[1].imshow(sedimemt, interpolation="none", cmap=colors.ListedColormap(['white', "gray"]))
        ax[1].set_title("sediment")
        img_dir = self._get_img_dir(parameters)
        fig.savefig(f"{img_dir}/poroisty and sediment volume contour.jpg", bbox_inches='tight', dpi=600)

        nx = eps.shape[0]
        ny = eps.shape[1]
        eps_internal_field = list()
        eps_internal_field.append(f"{str(nx*ny)}\n")
        eps_internal_field.append("(\n")
        for j in np.arange(0, ny):
            for i in np.arange(0, nx):
                eps_internal_field.append(f"{str(eps[i,j])}\n")
        eps_internal_field.append(")\n")
        eps_internal_field.append(";\n")

        # if not os.path.exists(f"{templateDir}"):
        #     os.makedirs(f"{templateDir}",exist_ok=True)

        templateDir = parameters["openfoam_file"]["openfoam_template_dir"]
        with open(f"{templateDir}/eps_template", "r") as fp:
            eps_template = fp.readlines()
        internal_field_line_index = 0
        for index, line in enumerate(eps_template):
            if line.startswith("internalField"):
                internal_field_line_index = index
                for i in range(len(eps_internal_field)):
                    eps_template.insert(internal_field_line_index+1+i, eps_internal_field[i])
                break

        targetDir = parameters["openfoam_file"]["openfoam_file_output_dir"]
        if not os.path.exists(f"{targetDir}/0"):
            os.makedirs(f"{targetDir}/0", exist_ok=True)
        if not os.path.exists(f"{targetDir}/system"):
            os.makedirs(f"{targetDir}/system", exist_ok=True)
        with open(f"{targetDir}/0/eps", "w") as fp:
            fp.writelines(eps_template)

        sediment_internal_field = list()
        sediment_internal_field.append(f"{str(nx*ny)}\n")
        sediment_internal_field.append("(\n")
        for j in np.arange(0, ny):
            for i in np.arange(0, nx):
                sediment_internal_field.append(f"{str(sedimemt[i,j])}\n")
        sediment_internal_field.append(")\n")
        sediment_internal_field.append(";\n")

        with open(f"{templateDir}/sediment_template", "r") as fp:
            sediment_template = fp.readlines()
        internal_field_line_index = 0
        for index, line in enumerate(sediment_template):
            if line.startswith("internalField"):
                internal_field_line_index = index
                for i in range(len(sediment_internal_field)):
                    sediment_template.insert(internal_field_line_index+1+i, sediment_internal_field[i])
                break

        with open(f"{targetDir}/0/sediment", "w") as fp:
            fp.writelines(sediment_template)

        with open(f"{templateDir}/blockMeshDict_template", "r") as fp:
            mesh_template = fp.readlines()
        for index, line in enumerate(mesh_template):
            if line.startswith("convertToMeters "):
                mesh_template.insert(index+1, f"lx {nx};\n")
                mesh_template.insert(index+2, f"ly {ny};\n")
                mesh_template.insert(index+3, f"Nx {nx};\n")
                mesh_template.insert(index+4, f"Ny {ny};\n")
                break

        with open(f"{targetDir}/system/blockMeshDict", "w") as fp:
            fp.writelines(mesh_template)
        print("OpenFOAM files were successfully outputted")

    def _select_random_deposition_subdomains(
            self, parameters, chunked_subdomain_nx, chunked_subdomain_ny, subdomains_dir):
        work_dir = parameters["work_dir"]
        chunked_subdomain_num_x = parameters["sediment_deposition"]["chunked_subdomain_num_x"]
        chunked_subdomain_num_y = parameters["sediment_deposition"]["chunked_subdomain_num_y"]
        min_grow_num = parameters["sediment_deposition"]["min_grow_num"]
        max_grow_num = parameters["sediment_deposition"]["max_grow_num"]
        seed = parameters["sediment_deposition"]["seed"]

        np.random.seed(seed)
        random_grow_numbers = np.random.randint(low=min_grow_num, high=max_grow_num+1,
                                                size=(chunked_subdomain_num_x, chunked_subdomain_num_y))  # user

        print(f"random_grow_number matrix are:\n {random_grow_numbers}")
        sediment_bearing_subdomains = []
        sediment_fraction = []
        for i in np.arange(0, chunked_subdomain_num_x):
            sediment_bearing_subdomains.append([])
            for j in np.arange(0, chunked_subdomain_num_y):
                k = chunked_subdomain_num_y*i+j
                subdomain_dir = os.path.join(subdomains_dir, f"subdomain{k}")
                os.chdir(subdomain_dir)
                grow_index = random_grow_numbers[i, j]
                if grow_index == 0:
                    sediment_bearing_subdomain = self._read_impl(
                        "./originalPoreSpace", chunked_subdomain_nx, chunked_subdomain_ny)
                    sediment_bearing_subdomain[sediment_bearing_subdomain == 1] = self._label_rock
                    sediment_bearing_subdomain[sediment_bearing_subdomain == 0] = self._label_pore
                else:
                    sediment_bearing_subdomain = self._read_synthetic(
                        "./originalPoreSpace", f"./randFilmTurnBack_grow_{grow_index}", chunked_subdomain_nx,
                        chunked_subdomain_ny)

                sediment_bearing_subdomains[i].append(sediment_bearing_subdomain)
                fraction = np.sum(sediment_bearing_subdomain ==
                                  self._label_sediment)/sediment_bearing_subdomain.size*100
                sediment_fraction.append(fraction)
                os.chdir(work_dir)

        print(f"mean sediment fraction: {round(np.mean(sediment_fraction),2)}%")
        print(f"std sediment fraction: {round(np.std(sediment_fraction),2)}%")

        fig, ax = plt.subplots()
        ax.bar(np.arange(0, len(sediment_fraction)), sediment_fraction)
        ax.set_ylabel("sediment volume fraction (%)")
        img_dir = self._get_img_dir(parameters)
        fig.savefig(f"{img_dir}/sediment_fraction_distribution.jpg", bbox_inches='tight', dpi=100)

        chunked_subdomain_num_x = parameters["sediment_deposition"]["chunked_subdomain_num_x"]
        chunked_subdomain_num_y = parameters["sediment_deposition"]["chunked_subdomain_num_y"]
        fig, axes = plt.subplots(
            chunked_subdomain_num_x, chunked_subdomain_num_y, figsize=(8, 8),
            subplot_kw={"xticks": [],
                        "yticks": []},
            gridspec_kw=dict(hspace=0.1, wspace=0.1))

        for i, ax in enumerate(axes.flat):
            ri = i//chunked_subdomain_num_y
            ci = i % chunked_subdomain_num_y
            sediment_bearing_subdomain = sediment_bearing_subdomains[ri][ci]
            grow_index = random_grow_numbers[ri, ci]
            if grow_index == 0:
                im = ax.imshow(sediment_bearing_subdomain, interpolation='none',
                               cmap=colors.ListedColormap(['black', 'gray']))
            else:
                im = ax.imshow(sediment_bearing_subdomain, interpolation='none',
                               cmap=colors.ListedColormap(['black', 'gray', 'white']))
        fig.savefig(f"{img_dir}/chunked_sediment_rock.jpg", bbox_inches='tight', dpi=100)
        return sediment_bearing_subdomains

    def _deposit_sediment_in_subdomains(
            self, parameters, chunked_subdomain_nx, chunked_subdomain_ny, subdomains_dir, subdomains):
        lsmpqs_bin_dir = parameters["sediment_deposition"]["lsmpqs_bin_dir"]
        max_grow_num = parameters["sediment_deposition"]["max_grow_num"]
        crystal_probability = parameters["sediment_deposition"]["crystal_probability"]
        work_dir = parameters["work_dir"]
        for i in np.arange(0, len(subdomains)):
            subdomain_dir = os.path.join(subdomains_dir, f"subdomain{i}")
            os.chdir(subdomain_dir)
            subprocess.check_output(
                f"{lsmpqs_bin_dir}/createMask2D subdomain{i}.raw {chunked_subdomain_ny} {chunked_subdomain_nx}",
                shell=True)
            subprocess.check_output(
                f"{lsmpqs_bin_dir}/synthesizeCokedStructure2D  1 {crystal_probability} {max_grow_num}", shell=True)
            os.chdir(work_dir)
        return max_grow_num, work_dir

    def _chunk_rock(self, parameters, raw_rock):
        domain_nx = parameters["rock_simulation"]["simulation_rock_nx"]
        domain_ny = parameters["rock_simulation"]["simulation_rock_ny"]
        chunked_subdomain_num_x = parameters["sediment_deposition"]["chunked_subdomain_num_x"]
        chunked_subdomain_num_y = parameters["sediment_deposition"]["chunked_subdomain_num_y"]

        if domain_nx % chunked_subdomain_num_x != 0:
            raise ValueError("chunked_subdomain_num_x should be evently divide by nx")

        if domain_ny % chunked_subdomain_num_y != 0:
            raise ValueError("chunked_subdomain_num_y should be evently divide by ny")

        chunked_subdomain_nx = int(domain_nx/chunked_subdomain_num_x)
        chunked_subdomain_ny = int(domain_ny/chunked_subdomain_num_y)

        work_dir = parameters["work_dir"]
        subdomains_dir = os.path.abspath(os.path.join(work_dir, "subdomains"))
        if not os.path.exists(subdomains_dir):
            os.mkdir(subdomains_dir)
        os.chdir(subdomains_dir)
        subdomains = self._chunk_impl(raw_rock, chunked_subdomain_nx, chunked_subdomain_ny)
        for i, subdomain in enumerate(subdomains):
            dir = os.path.abspath(f"./subdomain{i}")
            if not os.path.exists(dir):
                os.mkdir(dir)
            with open(os.path.join(dir, f"subdomain{i}.raw"), "w") as fp:
                subdomain.astype(np.uint8).tofile(fp)
        print(f"{len(subdomains)} subdomains were chunked")
        os.chdir(work_dir)

        fig, axes = plt.subplots(
            chunked_subdomain_num_x, chunked_subdomain_num_y, figsize=(8, 8),
            subplot_kw={"xticks": [],
                        "yticks": []},
            gridspec_kw=dict(hspace=0.1, wspace=0.1))

        for i, ax in enumerate(axes.flat):
            subdomain_dir = os.path.join(subdomains_dir, f"subdomain{i}")
            with open(os.path.join(subdomain_dir, f"subdomain{i}.raw"), "r") as fp:
                subdomain = np.fromfile(fp, dtype=np.uint8).reshape(chunked_subdomain_nx, chunked_subdomain_ny)
                im = ax.imshow(subdomain, interpolation='none', cmap=colors.ListedColormap(['black', 'gray']))

        img_dir = self._get_img_dir(parameters)
        fig.savefig(f"{img_dir}/chunked_raw_rock.jpg", bbox_inches='tight', dpi=100)
        return chunked_subdomain_nx, chunked_subdomain_ny, subdomains_dir, subdomains

    def _read_synthetic(
            self, origianl_fname, synthetic_fname, nx, ny, gb_size=3):
        original = self._read_impl(origianl_fname, nx, ny, gb_size)
        synthetic = self._read_impl(synthetic_fname, nx, ny, gb_size)
        sediment = synthetic-original
        composite = original*self._label_rock+sediment*self._label_sediment
        old_label_pore = 0
        composite[composite == old_label_pore] = self._label_pore
        return composite

    def _read_impl(self, fname, nx, ny, gb_size=3):
        with open(fname, "r") as fp:
            raw_data_gb = np.fromfile(fp, dtype=np.ubyte)
        raw_data_gb = raw_data_gb.reshape(nx+2*gb_size, ny+2*gb_size)
        raw_data = raw_data_gb[gb_size:nx+gb_size, gb_size:ny+gb_size]
        # raw data: solid is 0, void is 1
        # inverted raw data: solid is 1, void is 0
        raw_data_invert = np.zeros(raw_data.shape, dtype=np.uint8)
        raw_data_invert[raw_data == 0] = 1
        raw_data_invert[raw_data == 1] = 0
        return raw_data_invert

    def _show_rock(self, data, cmap=colors.ListedColormap(['black', 'gray', 'white']), ticks=[1, 2, 3], axis="off"):
        fig, ax = plt.subplots()
        im = ax.imshow(data, interpolation='none', cmap=cmap)
        ax.axis(axis)
        fig.colorbar(im, ticks=ticks)
        return ax, fig

    def _report_fraction(self, synthetic_data):
        print(f"raw rock porosity: {1-np.sum(synthetic_data==self._label_rock)/synthetic_data.size}")
        print(f"sediment-bearing rock porosity: {np.sum(synthetic_data==self._label_pore)/synthetic_data.size}")
        print(f"sediment volume fraction: {np.sum(synthetic_data==self._label_sediment)/synthetic_data.size}")

    def _chunk_impl(self, array, nrows, ncols):
        r, h = array.shape
        subdomains = []
        for subdomain in (array.reshape(h//nrows, nrows, -1, ncols)
                          .swapaxes(1, 2)
                          .reshape(-1, nrows, ncols)):
            subdomains.append(subdomain)
        return subdomains
