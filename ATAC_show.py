from tkinter import *
from tkinter.ttk import *
import tkinter as tk
import tkinter.ttk
from tkinter import filedialog
from PIL import Image, ImageTk
import os
from tkinter import Label
# from done import scRNA_sep
import tkinter.messagebox
import ATAC_sep


# 打开指定的图片文件，缩放至指定尺寸
def get_image(filename, width, height):
    im = Image.open(filename).resize((width, height))
    return ImageTk.PhotoImage(im)


def txt_read(files):
    txt_dict = {}
    fopen = open(files)
    for line in fopen.readlines():
        line = str(line).replace("\n", "")
        txt_dict[line.split(' ', 1)[0]] = line.split(' ', 1)[1]
        # split（）函数用法，逗号前面是以什么来分割，后面是分割成n+1个部分，且以数组形式从0开始
    fopen.close()
    return txt_dict


# txt_dict = txt_read('explain.txt')


class ATAC_App():
    def __init__(self, master):

        self.root = master
        self.root.title("scATAC-seq Task Progress")
        self.root.geometry("800x600")

        # 创建画布，设置要显示的图片，把画布添加至应用程序窗口
        canvas_root = tkinter.Canvas(self.root, width=800, height=600)
        im_root = get_image('窗口图片.jpg', 800, 600)
        canvas_root.create_image(400, 300, image=im_root)
        canvas_root.pack()

        # 初始化任务进度条
        self.progress = tk.DoubleVar()
        self.progress.set(0.0)

        # 创建label
        self.label = tk.Label(self.root, text="")
        # self.label.pack(padx=10, pady=10)
        self.label.place(x=30, y=10)

        global current_path
        current_path = os.getcwd()

        # 创建按钮并添加点击事件
        self.button_select = tk.Button(self.root, text="Select Folder", command=self.select_folder, fg="Gold",
                                       bg="#006400")  # 选择文件夹按钮
        # self.button_select.pack(padx=10, pady=10)
        self.button_select.place(x=30, y=40)

        global adata

        self.button1 = tk.Button(self.root, text="1.read data", command=self.do_task1, width=12)
        self.button1.place(x=100, y=120)
        self.button2 = tk.Button(self.root, text="2.TSS statistics", command=self.do_task2, width=14)
        self.button2.place(x=100, y=190)
        self.button3 = tk.Button(self.root, text="3.Quality control", command=self.do_task3, width=16)
        self.button3.place(x=100, y=260)
        self.button4 = tk.Button(self.root, text="4.Dimension reduction process", command=self.do_task4, width=28)
        self.button4.place(x=100, y=330)
        self.button5 = tk.Button(self.root, text="5.Find highly variable genes", command=self.do_task5, width=26)
        self.button5.place(x=320, y=120)
        self.button6 = tk.Button(self.root, text="6.Single gene plot", command=self.do_task6, width=18)
        self.button6.place(x=320, y=190)
        self.button7 = tk.Button(self.root, text="7.Cluster plot", command=self.do_task7, width=14)
        self.button7.place(x=320, y=260)

        self.show_select = tk.Button(self.root, text="Slide show result pictures", command=self.display_images,
                                     fg="OrangeRed", bg="Wheat")  # 幻灯片放映
        # self.show_select.pack(side = 'bottom')
        self.show_select.place(x=620, y=550)

        self.init()

        # 创建进度条并放置在窗口中
        # self.progressbar = Progressbar(self.root, orient="horizontal", length=200, mode="determinate", variable=self.progress)
        # self.progressbar.pack(pady=20)
        # self.progressbar.place(x=300,y=290)

        self.tip = tk.Label(self.root,
                            text="← Please select the data folder before pressing the data processing buttons.")
        self.tip.place(x=130, y=45)

        self.message = tk.Label(self.root,
                                text="Please follow the instructions.\nMake sure that the path you choose\nhas the following files:\n    filtered_feature_bc_matrix.h5\n    fragments.tsv.gz\n    singlecell.csv\n    .gtf", anchor=NW,justify='left')
        self.message.place(x=550, y=150)

        self.min_cells_message = tk.Label(self.root, text="read min cells:" ,anchor=NW,justify='right')
        self.min_cells_message.place(x=550, y=365)
        self.min_features_message = tk.Label(self.root, text="read min features:" ,anchor=NW,justify='right')
        self.min_features_message.place(x=550, y=410)

        self.min_cells = tk.Text(self.root, height=1.5, width=8)
        self.min_cells.place(x=680, y=360)
        self.min_features = tk.Text(self.root, height=1.5, width=8)
        self.min_features.place(x=680, y=410)

        self.back_main = tk.Button(self.root, text='back', font=('Helvetica', 15, 'bold'), command=self.back)
        self.back_main.place(x=20, y=550)

        self.root.mainloop()

    def init(self):
        self.button1.configure(state=DISABLED)
        self.button2.configure(state=DISABLED)
        self.button3.configure(state=DISABLED)
        self.button4.configure(state=DISABLED)
        self.button5.configure(state=DISABLED)
        self.button6.configure(state=DISABLED)
        self.button7.configure(state=DISABLED)

    def enable(self):
        self.button1.configure(state=NORMAL)
        self.button2.configure(state=NORMAL)
        self.button3.configure(state=NORMAL)
        self.button4.configure(state=NORMAL)
        self.button5.configure(state=NORMAL)
        self.button6.configure(state=NORMAL)
        self.button7.configure(state=NORMAL)

    def do_task1(self):

        self.button1.configure(text="busy...", state=DISABLED)
        if self.min_cells.get("1.0", "end") != '\n' :
            min_cells = int(self.min_cells.get("1.0", "end"))
        else:
            min_cells = 10
        if self.min_features.get("1.0", "end") != '\n':
            min_features = int(self.min_features.get("1.0", "end"))
        else:
            min_features = 200
        ATAC_sep.tran_1(folder_path,min_cells,min_features)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button1.configure(text="1.read data\n完成/重启", fg='DarkCyan', state=NORMAL)
        self.min_cells.destroy()
        self.min_features.destroy()
        self.min_cells_message.destroy()
        self.min_features_message.destroy()
        self.TSS_message = tk.Label(self.root, text="TSS.enrichment >", anchor=NW, justify='right')
        self.TSS_message.place(x=550, y=380)
        self.TSS = tk.Text(self.root, height=1.5, width=8)
        self.TSS.place(x=680, y=375)
        self.message.config(
            text="For scATAC-seq sequencing data, we recommend using the following QC indicators to evaluate the data quality:Nucleosome banding pattern,Transcriptional start site (TSS) enrichment score,Total number of fragments in peaks,Fraction of fragments in peaks",
            justify='left',wraplength = 200)
        self.message.update()

    def do_task2(self):

        self.button2.configure(text="busy...", state=DISABLED)
        if self.TSS.get("1.0", "end") != '\n':
            TSS_enrichment = int(self.TSS.get("1.0", "end"))
        else:
            TSS_enrichment = 1.2
        ATAC_sep.tran_2(folder_path,TSS_enrichment)
        self.TSS_message.destroy()
        self.TSS.destroy()
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button2.configure(text="2.TSS statistics\n完成/重启", fg='DarkCyan', state=NORMAL)


    def do_task3(self):

        self.button3.configure(text="busy...", state=DISABLED)
        ATAC_sep.tran_3(folder_path)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button3.configure(text="3.Quality control\n完成/重启", fg='DarkCyan', state=NORMAL)

        self.message.config(
            text="After linear dimensionality reduction, we can use the methods commonly used to analyze scRNA-seq data to cluster cells based on graph, and visualize the nonlinear dimensionality reduction (such as UMAP).",
            justify='left',wraplength = 200)

        self.dims_message = tk.Label(self.root, text="dims:", anchor=NW, justify='right')
        self.dims_message.place(x=590, y=365)
        self.k_param_message = tk.Label(self.root, text="k.param for FindNeighbors:", anchor=NW, justify='right')
        self.k_param_message.place(x=530, y=410)
        self.resolution_message = tk.Label(self.root, text="resolution for FindClusters:", anchor=NW, justify='right')
        self.resolution_message.place(x=530, y=455)


        self.dims = tk.Text(self.root, height=1.5, width=8)
        self.dims.place(x=680, y=360)
        self.k_param = tk.Text(self.root, height=1.5, width=8)
        self.k_param.place(x=680, y=410)
        self.resolution = tk.Text(self.root, height=1.5, width=8)
        self.resolution.place(x=680, y=460)

        self.message.update()

    def do_task4(self):

        self.button4.configure(text="busy...", state=DISABLED)

        # k.param=40,resolution=1,dims=30
        if self.dims.get("1.0", "end") != '\n':
            dims = int(self.dims.get("1.0", "end"))
        else:
            dims = 30
        if self.k_param.get("1.0", "end") != '\n':
            k_param = int(self.k_param.get("1.0", "end"))
        else:
            k_param = 40
        if self.resolution.get("1.0", "end") != '\n':
            resolution = int(self.resolution.get("1.0", "end"))
        else:
            resolution= 1

        ATAC_sep.tran_4(folder_path,k_param,resolution,dims)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button4.configure(text="4.Dimension reduction process\n完成/重启", fg='DarkCyan', state=NORMAL)

        self.message.config(
            text="Identify the 10 most highly variable genes.",
            justify='left',wraplength = 200)

        self.dims.destroy()
        self.dims_message.destroy()
        self.resolution.destroy()
        self.resolution_message.destroy()
        self.k_param.destroy()
        self.k_param_message.destroy()

        self.message.update()

    def do_task5(self):

        self.button5.configure(text="busy...", state=DISABLED)
        global top10_list
        top10,top10_list=ATAC_sep.tran_5(folder_path)

        self.root.update()
        tkinter.messagebox.showinfo('top10 gene id', top10)
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button5.configure(text="5.Find highly variable genes\n完成/重启", fg='DarkCyan', state=NORMAL)

        self.message.config(
            text="Check the differential expression of a single gene between different cluster cells and the expression distribution in all cells.",
            justify='left',wraplength = 200)

        self.gene_message = tk.Label(self.root, text="Which gene:", anchor=NW, justify='right')
        self.gene_message.place(x=550, y=365)

        self.gene = tk.Text(self.root, height=1.5, width=15)
        self.gene.place(x=680, y=360)

        self.message.update()

    def do_task6(self):
        global adata
        self.button6.configure(text="busy...", state=DISABLED)

        if self.gene.get("1.0", "end") != '\n':
            gene = self.gene.get("1.0", "end")
        else:
            gene=top10_list[1]

        ATAC_sep.tran_6(folder_path,gene)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button6.configure(text="6.Single gene plot\n完成/重启", fg='DarkCyan', state=NORMAL)

        self.message.config(
            text="Finally,view the distribution of all cell clusters.",
            justify='left')

        self.gene_message.destroy()
        self.gene.destroy()

        self.message.update()

    def do_task7(self):
        global adata
        self.button7.configure(text="busy...", state=DISABLED)
        ATAC_sep.tran_7(folder_path)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button7.configure(text="7.Cluster plot\n完成/重启", fg='DarkCyan', state=NORMAL)

    def select_folder(self):
        global folder_path
        folder_path = filedialog.askdirectory()
        self.label.config(text="Selected Folder: " + folder_path)
        self.tip.config(text="← data folder selected")
        # self.tip = tk.Label(self.root, text="← data folder selected")
        # self.tip.place(x=130, y=45)
        self.tip.update()
        self.enable()
        self.message.config(
            text="Now you can take the first step:\nread the data of the selected folder.\npress 1.read data",
            justify='left')
        self.message.update()

    def display_images(self):
        self.top = Toplevel()
        self.top.title('幻灯片展示')
        self.top.geometry("800x600")
        # 创建画布，设置要显示的图片，把画布添加至应用程序窗口
        canvas_top = tkinter.Canvas(self.top, width=800, height=600)
        im_top = get_image(current_path + '/窗口图片.jpg', 800, 600)
        canvas_top.create_image(400, 300, image=im_top)
        canvas_top.pack()
        images = []

        # upDir=os.path.abspath(os.path.join(os.path.dirname(folder_path), os.path.pardir))
        upDir = os.path.abspath(os.path.join(folder_path, "../"))
        print(upDir)
        folder_path_figures = upDir + r'\figure'
        print(folder_path_figures)
        for filename in os.listdir(folder_path_figures):
            if filename.endswith(('.jpg', '.png', '.jpeg')):
                image = Image.open(os.path.join(folder_path_figures, filename))
                image = image.resize((400, 400))
                images.append(ImageTk.PhotoImage(image))

        if images:
            index = 0

            def show_image():
                nonlocal index

                if index >= len(images):
                    # self.top.label.config(text="已经读取完毕")
                    # label(self.top, text="已经读取完毕")
                    widget = Label(self.top, text="已经读取完毕")
                    widget.config(bg='black', fg='yellow')

                    return
                # im_show = Label(self.top)
                # im_show.config(image=images[index])

                label_img = Label(self.top, image=images[index])
                label_img.configure(image=images[index])
                label_img.place(x=10, y=10)
                # label = tk.Label(self.top, image=images[index])
                # label.pack()

                widget = Label(self.top, text="This is the {} th picture.".format(index + 1), font=('Arial', 15))
                widget.config(bg='black', fg='yellow')
                widget.place(x=500, y=100)

                # explain = Label(self.top, text=txt_dict[str(index)], height=15, width=35,wraplength=300,font=('Arial', 12), justify='left')
                # explain.config(fg="#006400")
                # explain.place(x=450, y=150)

                # explain.update()

                index += 1

                self.top.after(4000, show_image)

            show_image()

        self.top.mainloop()

    def back(self):
        self.root.destroy()
        self.root = tk.Tk()
        from main_show import basedesk  # 防止互相调用进入死循环
        basedesk(self.root)

