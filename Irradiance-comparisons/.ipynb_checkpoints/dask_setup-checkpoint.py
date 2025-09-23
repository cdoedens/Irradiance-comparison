from dask.distributed import Client, wait
import psutil

def setup_dask_client(
    workload_type="io",        # "cpu", "io", or "mixed"
    max_workers=None,           # Limit max workers (default = all cores)
    reserve_mem_gb=50,          # Memory reserved for system and overhead
    max_mem_gb=None,            # Cap usable memory if needed
    dashboard=True
):
    """
    Setup Dask client with auto-scaling for workload type.

    Parameters:
        workload_type (str): Type of workload - "cpu", "io", or "mixed"
        max_workers (int): Max logical cores to use (default = system max)
        reserve_mem_gb (int): System memory to reserve
        max_mem_gb (int): Cap memory usage (defaults to system memory)
        dashboard (bool): Whether to print the dashboard link

    Returns:
        dask.distributed.Client
    """
    assert workload_type in ("cpu", "io", "mixed"), "Invalid workload_type"

    logical_cores = psutil.cpu_count(logical=True)
    total_memory_gb = psutil.virtual_memory().total / 1e9

    if max_workers is None:
        max_workers = logical_cores

    if max_mem_gb is None:
        max_mem_gb = total_memory_gb

    usable_mem_gb = max_mem_gb - reserve_mem_gb

    # Recommended presets based on workload type
    if workload_type == "cpu":
        threads_per_worker = 1
        n_workers = min(max_workers, logical_cores)
    elif workload_type == "io":
        threads_per_worker = 8
        n_workers = max(1, logical_cores // threads_per_worker)
    else:  # "mixed"
        threads_per_worker = 4
        n_workers = max(1, logical_cores // threads_per_worker)

    memory_per_worker = usable_mem_gb // n_workers

    client = Client(
        n_workers=n_workers,
        threads_per_worker=threads_per_worker,
        memory_limit=f"{int(memory_per_worker)}GB"
    )

    if dashboard:
        print(f"Dask dashboard: {client.dashboard_link}")

    return client