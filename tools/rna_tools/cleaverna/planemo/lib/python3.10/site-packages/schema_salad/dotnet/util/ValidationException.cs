namespace ${project_name};

public class ValidationException : Exception
{
    readonly string currentMessage;
    static readonly int indentPerLevel = 2;
    string bullet = "";
    private readonly List<ValidationException> children = new();

    public ValidationException() : this("")
    {
    }

    public ValidationException(string message)
        : this(message, new List<ValidationException>())
    {
    }

    public ValidationException(string message, ValidationException child)
        : this(message, new List<ValidationException> { child })
    {
    }

    public ValidationException(string message, List<ValidationException> children) : base(message)
    {
        currentMessage = message;
        foreach (ValidationException child in children)
        {
            this.children.AddRange(child.Simplify());
        }
    }

    private string Summary(int level)
    {
        string spaces = new(' ', level * indentPerLevel);
        return spaces + bullet + this.currentMessage;
    }

    public ValidationException WithBullet(string bullet)
    {
        this.bullet = bullet;
        return this;
    }

    List<ValidationException> Simplify()
    {
        if (this.ToString().Length > 0)
        {
            return new List<ValidationException> { this };
        }
        else return this.children;
    }

    private string PrettyString(int level = 0)
    {
        List<string> parts = new();
        int nextlevel = level;
        if (this.currentMessage.Length > 0)
        {
            parts.Add(this.Summary(level));
            nextlevel++;
        }

        foreach (ValidationException child in children)
        {
            parts.Add(child.PrettyString(nextlevel));
        }

        return string.Join('\n', parts);
    }

    public override string Message => PrettyString();

    public override string ToString()
    {
        return PrettyString();
    }
}
