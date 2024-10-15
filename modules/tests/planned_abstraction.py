# class Process:
#     def __init__(self, process_name, process_cmd, input_file=None, output_file=None, horizontal_name=None):
#         self.process_name = process_name
#         self.process_cmd = process_cmd
#         self.input_file = input_file
#         self.output_file = output_file
#         self.horizontal_name = horizontal_name
#
#     def run(self, self_app):
#         self_app.notify(f"{self.process_name} starting...", title=self.process_name)
#         logging.info(f"{self.process_name} starting...")
#         if self.horizontal_name:
#             self_app.query_one(self.horizontal_name).add_class("running")
#         threading.Thread(target=self._run, args=(self_app,)).start()
#
#     def _run(self, self_app):
#         try:
#             self._execute(self_app)
#             logging.info(f"{self.process_name} completed for {self.input_file}")
#             self_app.notify(f"{self.process_name} completed for {self.input_file}", title=self.process_name)
#         except Exception as e:
#             logging.error(f"An error occurred during {self.process_name}: {e}")
#             self_app.notify(
#                 f"An error occurred during {self.process_name}: {e}",
#                 severity="error",
#                 timeout=10.0,
#                 title=self.process_name,
#             )
#         finally:
#             if self.horizontal_name:
#                 self_app.query_one(self.horizontal_name).remove_class("running")
#
#     def _execute(self, self_app):
#         raise NotImplementedError("Subclasses must implement this method")
#
# class SamtoolsViewProcess(Process):
#     def __init__(self, input_file, output_file, view_region):
#         super().__init__(
#             "Samtools view",
#             f"samtools view {input_file} -o {output_file} region={view_region}",
#             input_file,
#             output_file,
#             "#sam_view_horizontal"
#         )
#         self.view_region = view_region
#
#     def _execute(self, self_app):
#         self_app.notify(self.process_cmd, title=self.process_name)
#         logging.info("Running command: " + self.process_cmd)
#
#         with pysam.AlignmentFile(self.input_file, "rb") as infile, pysam.AlignmentFile(
#             self.output_file, "wb", header=infile.header
#         ) as outfile:
#             for read in infile.fetch(region=self.view_region) if self.view_region else infile:
#                 outfile.write(read)
#
#
# def run_samtools_view(self):
#     if not self.workingDir:
#         self.workingDir = __file__
#     input_path = get_input(self, "sam_view_input_input")
#     output_path = get_input(self, "sam_view_output_input")
#     view_region = get_input(self, "sam_view_region_input")
#
#     if not input_path:
#         self.notify("Please provide a valid path", severity="warning", title="Samtools view")
#         return
#
#     if not output_path:
#         input_path_tuple = get_basename_and_ext(input_path)
#         last = len(input_path_tuple)-1
#         output_path = f"{input_path_tuple[0]}.view{input_path_tuple[last]}"
#
#     process = SamtoolsViewProcess(input_path, output_path, view_region)
#     process.run(self)
#
# class SamtoolsSortProcess(Process):
#     def __init__(self, input_file, output_file):
#         super().__init__(
#             "Samtools sort",
#             f"samtools sort {input_file} -o {output_file}",
#             input_file,
#             output_file,
#             "#sam_sort_horizontal"
#         )
#
#     def _execute(self, self_app):
#         self_app.notify(self.process_cmd, title=self.process_name)
#         logging.info("Running command: " + self.process_cmd)
#         pysam.sort("-o", self.output_file, self.input_file)
#
# class BcftoolsMpileupProcess(Process):
#     def __init__(self, input_file, reference, output_file):
#         super().__init__(
#             "Bcftools mpileup",
#             f"bcftools mpileup -f {reference} {input_file} -o {output_file}",
#             input_file,
#             output_file,
#             "#bcftools_mpileup_horizontal"
#         )
#         self.reference = reference
#
#     def _execute(self, self_app):
#         # Implement bcftools mpileup execution here
#         pass
#
# class VariantCallingPipeline:
#     def __init__(self, input_file, reference, output_prefix):
#         self.processes = [
#             SamtoolsSortProcess(input_file, f"{output_prefix}.sorted.bam"),
#             SamtoolsIndexProcess(f"{output_prefix}.sorted.bam"),
#             BcftoolsMpileupProcess(f"{output_prefix}.sorted.bam", reference, f"{output_prefix}.vcf"),
#             BcftoolsCallProcess(f"{output_prefix}.vcf", f"{output_prefix}.called.vcf"),
#         ]
#
#     def run(self, self_app):
#         for process in self.processes:
#             process.run(self_app)
#
# class Process:
#     def __init__(self, process_name, process_cmd, input_file=None, output_file=None, horizontal_name=None):
#         self.process_name = process_name
#         self.process_cmd = process_cmd
#         self.input_file = input_file
#         self.output_file = output_file
#         self.horizontal_name = horizontal_name
#
# def Process(self, process_name, process_cmd, input_id=None, output_id=None, horizontal_name=None):
#

