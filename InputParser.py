import os


class InputParser:
    def __init__(self, args, **kwargs):
        self.args = args
        self.path_to_R_main = os.path.join(os.path.dirname(os.path.realpath(__file__)), "from_immuneSIM/main.R")
        self.other_args = kwargs
        self.default_main_R_args = {
            "number_of_seqs": "number_of_sequences",
            "vdj_list": "list_germline_genes_allele_01",
            "species": "hs",
            "receptor": "ig",
            "chain": "h",
            "insertions_and_deletion_lengths": "insertions_and_deletion_lengths_df",
            "name_repertoire": "hs_igh_sim",
            "length_distribution_rand": "length_dist_simulation",
            "random": "FALSE",
            "shm_mode": 'none',
            "shm_prob": "15 / 350",
            "vdj_noise": "0",
            "vdj_dropout": "c(V=0, D=0, J=0)",
            "ins_del_dropout": "",
            "equal_cc": "FALSE",
            "freq_update_time": "round(0.5 * number_of_sequences)",
            "verbose": "TRUE",
            "airr_compliant": "TRUE"
        }
        self.need_raw_repr = {
            "number_of_seqs": False,
            "vdj_list": False,
            "species": True,
            "receptor": True,
            "chain": True,
            "insertions_and_deletion_lengths": False,
            "name_repertoire": True,
            "length_distribution_rand": False,
            "random": False,
            "shm.mode": True,
            "shm.prob": False,
            "vdj_noise": False,
            "vdj_dropout": False,
            "ins_del_dropout": True,
            "equal_cc": False,
            "freq_update_time": False,
            "verbose": False,
            "airr_compliant": False
        }

    def _check_args(self):
        assert (self.args.receptor == "ig" and self.args.chain in ("h", "l", "k")) or \
               (self.args.receptor == "tr" and self.args.chain in ("a", "b")), \
            f"Check receptor and chain pair because the one given ({self.args.receptor} and {self.args.chain}) can not be processed"

    def _parse_python_input_to_mainR(self):
        """kwargs are meant for manually modifying non-essential input field for simulation.
        'smh.mode' and 'smh.prob' should be passed as 'smh_mode' and 'smh_prob'
        Ex:
            freq_update_time = "round(0.5 * number_of_sequences)"
            smh_mode = "naive"
            shm.prob = "15/350"
        """
        essential_mainR_field = {
            "species": f"{self.args.species}",
            "receptor": f"{self.args.receptor}",
            "chain": f"{self.args.chain}",
            "name_repertoire": f"{self.args.name_repertoire}" if self.args.name_repertoire != "NA"
            else f"{self.args.species}_{self.args.receptor}{self.args.chain}_sim"
        }
        self._check_args()

        for parameter, value in essential_mainR_field.items():
            self.default_main_R_args[parameter] = value

        for parameter, value in self.other_args.items():
            self.default_main_R_args[parameter] = value

        # need to rename shm_mode and shm_prob to shm.mode shm.prob
        self.default_main_R_args['shm.mode'] = self.default_main_R_args.pop('shm_mode')
        self.default_main_R_args['shm.prob'] = self.default_main_R_args.pop('shm_prob')

    def _modify_mainR_with_params(self):
        with open(self.path_to_R_main, "r") as cfile:
            raw_mainR_content: str = cfile.read()
        split_by_lines_mainR_content: list[str] = raw_mainR_content.split("\n")

        # ======== add the given inputs to main.R =========
        split_by_lines_mainR_content[1] = f"number_of_sequences <- {self.args.number_of_seqs}"
        default_args_dump_container = [
            f"{' ' * 21}{key} = {value if not self.need_raw_repr[key] else repr(value)}"
            f"{')' if idx == len(self.default_main_R_args) - 1 else ','}"
            for idx, (key, value) in enumerate(self.default_main_R_args.items())]
        split_by_lines_mainR_content[3:3 + len(default_args_dump_container)] = default_args_dump_container
        with open(self.path_to_R_main, "w") as cfile:
            cfile.write("\n".join(split_by_lines_mainR_content))

    def parse_and_modify_mainR(self):
        self._parse_python_input_to_mainR()
        self._modify_mainR_with_params()
