import sys
import os


def main():
    IPYNB_FILENAME = 'render_pdb.ipynb'
    CONFIG_FILENAME = '.config_ipynb'
    with open(CONFIG_FILENAME, 'w') as f:
        f.write(' '.join(sys.argv))
    os.system('jupyter nbconvert --execute {:s} --to html'.format(IPYNB_FILENAME))
    os.system(IPYNB_FILENAME.replace(".ipynb", ".html"))


if __name__ == "__main__":
    main()