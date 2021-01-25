import matplotlib.pyplot as plt
import seaborn as sns


def set_style_paper():
    sns.set_context('paper')
    sns.set_style(style='ticks')
    fontsize = 7
    plt.rcParams['figure.figsize'] = [2.67, 2]
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['legend.fontsize'] = fontsize
    plt.rcParams['legend.title_fontsize'] = fontsize
    plt.rcParams['figure.titlesize'] = fontsize
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['svg.fonttype'] = 'none'


def set_style_talk():
    sns.set_context('talk')
    sns.set_style(style='ticks')
    plt.rcParams['svg.fonttype'] = 'none'


# palettes
saturation = 0.75
colors = 8
palettes = dict(T7=sns.dark_palette(sns.desaturate(color=(1, 0.5, 0), prop=saturation),
                                    reverse=True, as_cmap=False, n_colors=colors, input='rgb'),
                ORC=sns.dark_palette(sns.desaturate(color=(0, 1, 0.5), prop=saturation),
                                     reverse=True, as_cmap=False, n_colors=colors, input='rgb'),
                MCM=sns.dark_palette(sns.desaturate(color=(0, 0.5, 1), prop=saturation),
                                     reverse=True, as_cmap=False, n_colors=colors, input='rgb'),
                NUC=sns.dark_palette(sns.desaturate(color=(1, 0, 0.5), prop=saturation),
                                     reverse=True, as_cmap=False, n_colors=colors, input='rgb'),
                MCM_NUC=sns.dark_palette(sns.desaturate(color=(0.5, 0.25, 0.75), prop=saturation),
                                         reverse=True, as_cmap=False, n_colors=8, input='rgb'),
                OCCM=sns.dark_palette(sns.desaturate(color=(0, 0.75, 0.75), prop=saturation),
                                      reverse=True, as_cmap=False, n_colors=8, input='rgb'),
                qualitative=sns.color_palette('Set2', 8))
