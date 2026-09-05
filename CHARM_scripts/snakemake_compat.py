"""Enable recorded parameter/input invalidation in the pinned Snakemake 5.20."""

from snakemake.jobs import Job
from snakemake.rules import Rule


def track_rule_changes(workflow):
    # Input functions run while building the DAG, before update_needrun().
    # Reuse native persistence comparisons without writing state during dry-runs.
    def change_check(rule):
        def check(wildcards):
            persistence = workflow.persistence
            if persistence is None:
                return []
            job = Job(rule, persistence.dag, wildcards_dict=dict(wildcards.items()))
            for output in job.output:
                # DAG construction also tries wildcard matches that lose to a
                # more specific producer (e.g. count versus convertCountFormat).
                if persistence.metadata(output).get("rule") != rule.name:
                    continue
                if (
                    persistence.params_changed(job, file=output)
                    or persistence.input_changed(job, file=output)
                    or persistence.code_changed(job, file=output)
                ):
                    persistence.dag.forcefiles.add(str(output))
                    break
            return []
        return check

    for rule in workflow.rules:
        if rule.output:
            rule.set_input(change_check(Rule(rule)))
