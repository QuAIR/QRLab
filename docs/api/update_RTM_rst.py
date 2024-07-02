"""
How to use?
Run this python file on the root directory of this repository (Attention: Not its subdirectory) and then all files in source will be generated.
"""

import os

# 项目名称
_platform_name = "QRLab"

# 忽略的文件列表
_ignore_file_names = [
    ".git",
    "_",
    "docs",
    "entanglement_theory.Fawzi",
    "entanglement_theory.SeesawLOCC",
    "static_coherence.output_directory",
    "static_coherence.testing",
    "Magic.Magic_Qubit.Test_RobMag",
    "test_files",
]

# Sphinx文档源目录
_sphinx_source_dir = os.path.join(".", "docs", "api", "sphinx_src")


def is_correct_directory():
    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)
    required_file = os.path.join(script_dir, "docs", "api")
    return os.path.exists(required_file)


def _list_matlab_files(path=".", base_path="", file_name_attr_list=None):
    """
    递归列出给定目录中的文件和文件夹。

    Args:
        path (str): 开始列出文件的目录路径。
        base_path (str, optional): 相对路径计算的基路径。默认为空字符串。
        file_name_attr_list (list, optional): 用于存储结果的列表。默认为None。

    Returns:
        list: 包含每个文件的相对路径和类型（文件夹或MATLAB文件）的元组列表。
    """
    if file_name_attr_list is None:
        file_name_attr_list = []

    for child in os.listdir(path):
        if child.startswith("__"):
            continue
        child_path = os.path.join(path, child)
        relative_path = os.path.join(base_path, child).replace(os.path.sep, ".")

        if os.path.isdir(child_path):
            # 递归获取子目录内容
            sub_list = _list_matlab_files(child_path, relative_path, [])
            if sub_list:  # 仅当子目录非空时才添加
                file_name_attr_list.append((relative_path, "folder"))
                file_name_attr_list.extend(sub_list)
        elif child.endswith(".m"):
            file_name_attr_list.append((relative_path[:-2], "file"))  # 去掉.m扩展名
    file_name_attr_list = [
        sub_array
        for sub_array in file_name_attr_list
        if not any(
            sub_array[0].startswith(ignore_item) for ignore_item in _ignore_file_names
        )
    ]
    file_name_attr_list.sort()
    return file_name_attr_list


def _update_index_rst(file_name_attr_list):
    """
    更新index.rst文件。

    Args:
        file_name_attr_list (list): 文件名属性列表。
    """
    title = f":hide-footer:\nWelcome to {_platform_name}'s documentation!"
    content = f"{title}\n{'=' * len(title)}\n\n"
    content += """\
`Go to QuAIR Home <https://quair.github.io/>`_

.. automodule:: .

.. toctree::
   :maxdepth: 1

"""

    folders = [
        name
        for name, attr in file_name_attr_list
        if attr == "folder" and "." not in name
    ]
    for folder in folders:
        content += f"   {folder}\n"

    with open(os.path.join(_sphinx_source_dir, "index.rst"), "w") as file:
        file.write(content)


def _update_folder_rst(folder_name, file_name_attr_list):
    """
    为文件夹生成rst文件。

    Args:
        folder_name (str): 文件夹名称。
        file_name_attr_list (list): 文件名属性列表。
    """
    content = f"{folder_name}\n{'=' * len(folder_name)}\n\n.. automodule:: {folder_name}\n\n.. toctree::\n   :maxdepth: 1\n\n"

    # 获取子文件夹
    sub_folders = [
        name
        for name, attr in file_name_attr_list
        if attr == "folder"
        and name.startswith(folder_name)
        and name.count(".") == folder_name.count(".") + 1
    ]
    # 获取当前目录下的MATLAB文件
    matlab_files = [
        name
        for name, attr in file_name_attr_list
        if attr == "file"
        and name.startswith(folder_name)
        and name.count(".") == folder_name.count(".") + 1
    ]

    for sub_folder in sub_folders:
        content += f"   {sub_folder}\n"

    for matlab_file in matlab_files:
        content += f"\n.. autofunction:: {matlab_file}\n"

    with open(os.path.join(_sphinx_source_dir, f"{folder_name}.rst"), "w") as file:
        file.write(content)


def _update_conf_py():
    """
    更新Sphinx配置文件conf.py。
    """
    conf_content = """
import os
import sys

sys.path.insert(0, os.path.join("..", "..", ".."))
matlab_src_dir = os.path.abspath(os.path.join("..", "..", ".."))
extensions = [
    "sphinx.ext.autodoc",
    "sphinxcontrib.matlab",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx_immaterial",
]
primary_domain = "mat"
# project = "QRLab"
master_doc = "index"
source_suffix = ".rst"
nitpicky = True


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = "sphinx_rtd_theme"
html_theme = "sphinx_immaterial"
html_title = "QRLab"
html_short_title = "QRLab"
build_dir = "api"
# html_theme_options = {
#     'navigation_depth': 1,
# }
html_theme_options = {
    "base_url": "https://quair.github.io/QRLab/",
    "repo_url": "https://github.com/QuAIR/QRLab",
    "repo_name": "QRLab",
    # 'google_analytics_account': 'UA-XXXXX',
    "html_minify": True,
    "css_minify": True,
    "nav_title": "QRLab API Documentation",
    # 'logo_icon': '&#xe869',
    # 'globaltoc_depth': 2,
    "palette": { "primary": "orange" }
}
html_favicon = '../favicon.svg'
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

master_doc = "index"

# Autodoc
napoleon_numpy_docstring = False
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_warningiserror = False
autodoc_inherit_docstrings = False
autodoc_docstring_signature = False
autodoc_typehints_description_target = "documented"
autodoc_typehints_format = "short"
"""

    file_path = os.path.join(_sphinx_source_dir, "conf.py")
    with open(file_path, "w") as file:
        file.write(conf_content)


if __name__ == "__main__":
    _current_script_path = os.path.abspath(__file__)
    _platform_dir_path = os.path.dirname(
        os.path.dirname(os.path.dirname(_current_script_path))
    )
    _current_working_dir = os.getcwd()
    if _current_working_dir == _platform_dir_path:
        result = _list_matlab_files()
        os.makedirs(_sphinx_source_dir, exist_ok=True)
        _update_index_rst(result)
        for name, attr in result:
            if attr == "folder":
                _update_folder_rst(name, result)
        _update_conf_py()
    else:
        raise SystemExit(f"The current working directory is not {_platform_dir_path}.")
